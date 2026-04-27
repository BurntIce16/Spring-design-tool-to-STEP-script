#!/usr/bin/env python3
'''
build_spring_cam_v4.py

Single-file helper for spring-design-tool CSV output.

Compared with v3, this version can also export STEP using B-spline curves and simplifies the CSV polylines before creating
OpenSCAD/STEP geometry. This produces much smaller and more usable STEP files, and can avoid one STEP edge per CSV segment.

It:
  1. Reads outer_m.csv, inner_m.csv, and cam profile CSV files.
  2. Converts coordinates from meters to millimeters.
  3. Auto-repairs local self-intersecting loops in cam_profile_m.csv when needed.
  4. Simplifies closed polylines using Ramer-Douglas-Peucker tolerance in mm.
  5. Writes simplified CSVs into the output directory.
  6. Generates an OpenSCAD file from the simplified profiles.
  7. Optionally opens the generated .scad file in OpenSCAD.
  8. Optionally exports real STEP files using FreeCADCmd / FreeCAD.

Default expected files in the current directory:
  outer_m.csv
  inner_m.csv
  cam_profile_repaired_m.csv  preferred if present
  cam_profile_m.csv           fallback; auto-repaired by default

Typical usage:
  python3 build_spring_cam_v4.py --separate-step

More conservative:
  python3 build_spring_cam_v4.py --separate-step --simplify-tolerance 0.005

More aggressive:
  python3 build_spring_cam_v4.py --separate-step --simplify-tolerance 0.02

Disable simplification:
  python3 build_spring_cam_v4.py --separate-step --simplify-tolerance 0

Dependencies:
  - OpenSCAD for opening/rendering the generated .scad file.
  - FreeCADCmd or freecadcmd for true STEP export.
'''

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable


Point = tuple[float, float]


def read_points_csv(path: Path, scale: float = 1000.0) -> list[Point]:
    points: list[Point] = []

    with path.open(newline='') as f:
        reader = csv.reader(f)
        for row_number, row in enumerate(reader, start=1):
            if not row:
                continue

            try:
                x = float(row[0]) * scale
                y = float(row[1]) * scale
            except (IndexError, ValueError) as e:
                raise ValueError(f'{path}:{row_number}: expected at least two numeric CSV columns') from e

            points.append((x, y))

    if len(points) < 3:
        raise ValueError(f'{path}: expected at least 3 points, got {len(points)}')

    return points


def write_points_csv_mm_as_m(path: Path, points_mm: list[Point]) -> None:
    # The spring-design-tool CSV convention is meters.
    with path.open('w', newline='') as f:
        writer = csv.writer(f)
        for x_mm, y_mm in points_mm:
            writer.writerow([f'{x_mm / 1000.0:.18g}', f'{y_mm / 1000.0:.18g}'])


def remove_duplicate_closure(points: list[Point], eps: float = 1e-9) -> tuple[list[Point], bool]:
    if len(points) > 1 and math.dist(points[0], points[-1]) <= eps:
        return points[:-1], True
    return list(points), False


def close_points(points: list[Point]) -> list[Point]:
    if not points:
        return points

    if math.dist(points[0], points[-1]) <= 1e-9:
        return points

    return points + [points[0]]


def point_segment_distance(p: Point, a: Point, b: Point) -> float:
    ax, ay = a
    bx, by = b
    px, py = p

    dx = bx - ax
    dy = by - ay
    length2 = dx * dx + dy * dy

    if length2 <= 0.0:
        return math.dist(p, a)

    t = ((px - ax) * dx + (py - ay) * dy) / length2

    if t <= 0.0:
        return math.dist(p, a)

    if t >= 1.0:
        return math.dist(p, b)

    qx = ax + t * dx
    qy = ay + t * dy
    return math.hypot(px - qx, py - qy)


def rdp_open(points: list[Point], tolerance_mm: float) -> list[Point]:
    if tolerance_mm <= 0.0 or len(points) <= 2:
        return list(points)

    keep = [False] * len(points)
    keep[0] = True
    keep[-1] = True

    stack = [(0, len(points) - 1)]

    while stack:
        start, end = stack.pop()
        a = points[start]
        b = points[end]

        max_distance = -1.0
        max_index = -1

        for i in range(start + 1, end):
            distance = point_segment_distance(points[i], a, b)
            if distance > max_distance:
                max_distance = distance
                max_index = i

        if max_distance > tolerance_mm:
            keep[max_index] = True
            stack.append((start, max_index))
            stack.append((max_index, end))

    return [point for point, should_keep in zip(points, keep) if should_keep]


def simplify_closed_polyline(points: list[Point], tolerance_mm: float) -> list[Point]:
    '''
    Simplifies a closed polyline without treating closure as a single degenerate
    segment. It splits the loop into two open chains and simplifies each chain.
    '''
    open_points, _was_closed = remove_duplicate_closure(points)

    if tolerance_mm <= 0.0 or len(open_points) <= 3:
        return close_points(open_points)

    # Split at approximately opposite points.
    p0 = open_points[0]
    split_index = max(range(1, len(open_points)), key=lambda i: math.dist(p0, open_points[i]))

    chain_a = open_points[: split_index + 1]
    chain_b = open_points[split_index:] + [open_points[0]]

    simplified_a = rdp_open(chain_a, tolerance_mm)
    simplified_b = rdp_open(chain_b, tolerance_mm)

    combined = simplified_a[:-1] + simplified_b[:-1]

    # Remove accidental consecutive duplicates.
    cleaned: list[Point] = []
    for point in combined:
        if not cleaned or math.dist(cleaned[-1], point) > 1e-9:
            cleaned.append(point)

    if len(cleaned) < 3:
        return close_points(open_points)

    return close_points(cleaned)


def segment_intersection(
    a: Point,
    b: Point,
    c: Point,
    d: Point,
    eps: float = 1e-12,
) -> Point | None:
    r = (b[0] - a[0], b[1] - a[1])
    s = (d[0] - c[0], d[1] - c[1])

    denom = r[0] * s[1] - r[1] * s[0]
    if abs(denom) <= eps:
        return None

    qmp = (c[0] - a[0], c[1] - a[1])
    t = (qmp[0] * s[1] - qmp[1] * s[0]) / denom
    u = (qmp[0] * r[1] - qmp[1] * r[0]) / denom

    if eps < t < 1.0 - eps and eps < u < 1.0 - eps:
        return (a[0] + t * r[0], a[1] + t * r[1])

    return None


def find_first_self_intersection(points: list[Point]) -> tuple[int, int, Point] | None:
    open_points, _closed = remove_duplicate_closure(points)
    n = len(open_points)

    for i in range(n):
        a = open_points[i]
        b = open_points[(i + 1) % n]

        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue

            c = open_points[j]
            d = open_points[(j + 1) % n]

            intersection = segment_intersection(a, b, c, d)
            if intersection is not None:
                return i, j, intersection

    return None


def count_self_intersections(points: list[Point], max_count: int | None = None) -> int:
    open_points, _closed = remove_duplicate_closure(points)
    n = len(open_points)
    count = 0

    for i in range(n):
        a = open_points[i]
        b = open_points[(i + 1) % n]

        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue

            c = open_points[j]
            d = open_points[(j + 1) % n]

            if segment_intersection(a, b, c, d) is not None:
                count += 1
                if max_count is not None and count >= max_count:
                    return count

    return count


def repair_local_self_intersections(points: list[Point]) -> tuple[list[Point], int]:
    '''
    Removes small loop artifacts from a closed polyline.

    For an intersection between segment i->i+1 and segment j->j+1, the loop
    between i+1 and j is removed and replaced by the intersection point.

    This is intentionally conservative and intended for the spring-design-tool
    cam_profile_m.csv artifacts observed here. It is not a general-purpose CAD
    healing system.
    '''
    repaired, _closed = remove_duplicate_closure(points)
    repairs = 0

    while True:
        hit = find_first_self_intersection(repaired)
        if hit is None:
            break

        i, j, intersection = hit
        repaired = repaired[: i + 1] + [intersection] + repaired[j + 1 :]
        repairs += 1

        if repairs > 1000:
            raise RuntimeError('too many cam repair iterations; aborting')

    return close_points(repaired), repairs


def simplify_and_report(name: str, points: list[Point], tolerance_mm: float) -> list[Point]:
    before = len(points)
    simplified = simplify_closed_polyline(points, tolerance_mm)
    after = len(simplified)

    intersections = count_self_intersections(simplified, max_count=1)
    if intersections:
        raise RuntimeError(f'{name}: simplification created or preserved a self-intersection')

    if tolerance_mm > 0.0:
        print(f'{name}: simplified {before} -> {after} points at tolerance {tolerance_mm:g} mm')
    else:
        print(f'{name}: simplification disabled, {after} points')

    return simplified


def find_executable(explicit: str | None, candidates: Iterable[str]) -> str | None:
    if explicit:
        return explicit

    for candidate in candidates:
        found = shutil.which(candidate)
        if found:
            return found

    return None


def format_scad_points(name: str, points: list[Point]) -> str:
    lines = [f'{name} = [']
    for x, y in points:
        lines.append(f'  [{x:.9f}, {y:.9f}],')
    lines.append('];')
    return '\n'.join(lines)


def write_scad(
    path: Path,
    outer_points: list[Point],
    inner_points: list[Point],
    cam_points: list[Point],
    spring_height_mm: float,
    cam_height_mm: float,
) -> None:
    scad = f'''// Generated by build_spring_cam_v4.py
//
// Bodies:
//   spring body = simplified outer profile - simplified inner profile
//   cam body    = simplified selected/repaired cam profile
//
// Coordinates are converted from meters to millimeters by the Python script.

$fn = 48;

spring_height_mm = {spring_height_mm:.9f};
cam_height_mm = {cam_height_mm:.9f};

// "assembly", "spring", "cam", "outlines"
export_mode = "assembly";

// Z offset is only for visual inspection.
// Leave 0.0 for coaxial placement.
cam_z_offset_mm = 0.0;

// Usually leave this at zero because Python already simplified the profiles.
profile_cleanup_mm = 0.000;

outline_width_mm = 0.035;

module poly(points) {{
    if (profile_cleanup_mm == 0)
        polygon(points = points);
    else
        offset(r = profile_cleanup_mm)
            offset(r = -profile_cleanup_mm)
                polygon(points = points);
}}

module spring_body() {{
    linear_extrude(height = spring_height_mm)
        difference() {{
            poly(outer_points);
            poly(inner_points);
        }}
}}

module cam_body() {{
    translate([0, 0, cam_z_offset_mm])
        linear_extrude(height = cam_height_mm)
            poly(cam_points);
}}

module outline_curve(points, w = outline_width_mm, closed = true) {{
    n = len(points);

    for (i = [0 : n - 2])
        hull() {{
            translate(points[i]) circle(r = w);
            translate(points[i + 1]) circle(r = w);
        }}

    if (closed && n > 2)
        hull() {{
            translate(points[n - 1]) circle(r = w);
            translate(points[0]) circle(r = w);
        }}
}}

{format_scad_points("outer_points", outer_points)}

{format_scad_points("inner_points", inner_points)}

{format_scad_points("cam_points", cam_points)}

echo(str("outer point count = ", len(outer_points)));
echo(str("inner point count = ", len(inner_points)));
echo(str("cam point count = ", len(cam_points)));
echo(str("export_mode = ", export_mode));

if (export_mode == "assembly") {{
    color("lightgray") render(convexity = 10) spring_body();
    color("orange") render(convexity = 10) cam_body();
}}
else if (export_mode == "spring")
    render(convexity = 10) spring_body();
else if (export_mode == "cam")
    render(convexity = 10) cam_body();
else if (export_mode == "outlines") {{
    outline_curve(outer_points);
    outline_curve(inner_points);
    outline_curve(cam_points);
}}
else
    echo("ERROR: export_mode must be assembly, spring, cam, or outlines");
'''
    path.write_text(scad)


def open_scad_gui(openscad: str, scad_path: Path) -> None:
    subprocess.Popen([openscad, str(scad_path)])


def export_mesh_with_openscad(
    openscad: str,
    scad_path: Path,
    output_path: Path,
    export_mode: str,
) -> None:
    subprocess.run(
        [
            openscad,
            '-o',
            str(output_path),
            '-D',
            f'export_mode="{export_mode}"',
            str(scad_path),
        ],
        check=True,
    )


def write_freecad_export_script(
    path: Path,
    outer_csv: Path,
    inner_csv: Path,
    cam_csv: Path,
    output_step: Path,
    spring_output_step: Path | None,
    cam_output_step: Path | None,
    spring_height_mm: float,
    cam_height_mm: float,
    step_curve_mode: str,
    spline_approx_tolerance_mm: float,
    spline_break_angle_deg: float,
) -> None:
    # This script is executed by FreeCADCmd, not by normal Python.
    freecad_script = f'''import csv
import math
from pathlib import Path

import FreeCAD
import Part
import Import


SCALE = 1000.0

OUTER_CSV = Path(r"{outer_csv}")
INNER_CSV = Path(r"{inner_csv}")
CAM_CSV = Path(r"{cam_csv}")

OUTPUT_STEP = Path(r"{output_step}")
SPRING_OUTPUT_STEP = {repr(str(spring_output_step)) if spring_output_step else "None"}
CAM_OUTPUT_STEP = {repr(str(cam_output_step)) if cam_output_step else "None"}

SPRING_HEIGHT_MM = {spring_height_mm!r}
CAM_HEIGHT_MM = {cam_height_mm!r}

STEP_CURVE_MODE = "{step_curve_mode}"
SPLINE_APPROX_TOLERANCE_MM = {spline_approx_tolerance_mm!r}
SPLINE_BREAK_ANGLE_DEG = {spline_break_angle_deg!r}


def read_points(path):
    points = []
    with path.open(newline="") as f:
        for row_number, row in enumerate(csv.reader(f), start=1):
            if not row:
                continue
            try:
                points.append(FreeCAD.Vector(float(row[0]) * SCALE, float(row[1]) * SCALE, 0))
            except (IndexError, ValueError) as e:
                raise RuntimeError(f"{{path}}:{{row_number}}: expected at least two numeric columns") from e

    if len(points) < 3:
        raise RuntimeError(f"{{path}}: expected at least 3 points, got {{len(points)}}")

    if points[0].distanceToPoint(points[-1]) > 1e-7:
        points.append(points[0])

    return points


def strip_duplicate_closure(points):
    if len(points) > 1 and points[0].distanceToPoint(points[-1]) <= 1e-7:
        return points[:-1]
    return list(points)


def make_polyline_wire(points):
    return Part.Wire(Part.makePolygon(points))


def make_line_edge(a, b):
    return Part.makeLine(a, b)


def edge_points_are_closed(edge):
    try:
        return edge.Vertexes[0].Point.distanceToPoint(edge.Vertexes[-1].Point) <= 1e-7
    except Exception:
        return False


def close_edges_if_needed(edges):
    wire = Part.Wire(edges)
    if wire.isClosed():
        return wire

    try:
        start = edges[0].Vertexes[0].Point
        end = edges[-1].Vertexes[-1].Point
    except Exception as e:
        raise RuntimeError("failed to determine wire endpoints for closure") from e

    if start.distanceToPoint(end) > 1e-7:
        edges.append(make_line_edge(end, start))

    return Part.Wire(edges)


def make_open_spline_edge(name, points):
    # For very short spans a spline is not worth it and may fail.
    if len(points) < 3:
        return make_line_edge(points[0], points[-1])

    last_error = None

    if STEP_CURVE_MODE == "spline-approx":
        curve = Part.BSplineCurve()
        try:
            curve.approximate(
                Points=points,
                Tolerance=SPLINE_APPROX_TOLERANCE_MM,
                DegMin=3,
                DegMax=5,
            )
            return curve.toShape()
        except Exception as e:
            last_error = e
            # Fall through to interpolation.

    curve = Part.BSplineCurve()
    try:
        curve.interpolate(points)
        return curve.toShape()
    except Exception as e:
        last_error = e

    # Final fallback: use a small polyline for this span.
    try:
        return Part.Wire(Part.makePolygon(points))
    except Exception as e:
        raise RuntimeError(f"{{name}}: failed to create spline or fallback polyline span") from last_error or e


def make_closed_periodic_spline_wire(name, points):
    pts = strip_duplicate_closure(points)
    last_error = None

    # For a fully smooth loop, use one periodic B-spline edge.
    if len(pts) >= 4:
        if STEP_CURVE_MODE == "spline-approx":
            curve = Part.BSplineCurve()
            try:
                curve.approximate(
                    Points=pts + [pts[0]],
                    Tolerance=SPLINE_APPROX_TOLERANCE_MM,
                    DegMin=3,
                    DegMax=5,
                )
                edge = curve.toShape()
                wire = Part.Wire([edge])
                if wire.isClosed():
                    return wire
            except Exception as e:
                last_error = e

        # Interpolation supports PeriodicFlag on FreeCAD builds using OCC.
        attempts = [
            lambda c: c.interpolate(Points=pts, PeriodicFlag=True),
            lambda c: c.interpolate(pts, True),
            lambda c: c.interpolate(pts + [pts[0]]),
        ]

        for attempt in attempts:
            curve = Part.BSplineCurve()
            try:
                attempt(curve)
                edge = curve.toShape()
                wire = Part.Wire([edge])
                if wire.isClosed():
                    return wire

                # Non-periodic fallback with explicit closure.
                return close_edges_if_needed([edge])
            except Exception as e:
                last_error = e

    raise RuntimeError(f"{{name}}: failed to create closed spline wire") from last_error


def angle_between(v1, v2):
    n1 = v1.Length
    n2 = v2.Length

    if n1 <= 1e-12 or n2 <= 1e-12:
        return 0.0

    dot = max(-1.0, min(1.0, v1.dot(v2) / (n1 * n2)))
    return math.degrees(math.acos(dot))


def find_break_indices(points):
    pts = strip_duplicate_closure(points)
    breaks = []
    n = len(pts)

    for i in range(n):
        prev_point = pts[(i - 1) % n]
        point = pts[i]
        next_point = pts[(i + 1) % n]

        v1 = point.sub(prev_point)
        v2 = next_point.sub(point)

        turn = angle_between(v1, v2)
        if turn >= SPLINE_BREAK_ANGLE_DEG:
            breaks.append(i)

    # Avoid generating hundreds of tiny spline spans if the threshold is too low.
    if len(breaks) > max(16, n // 8):
        return []

    return breaks


def make_segmented_spline_wire(name, points):
    pts = strip_duplicate_closure(points)
    breaks = find_break_indices(points)

    if not breaks:
        return make_closed_periodic_spline_wire(name, points)

    edges = []
    break_count = len(breaks)

    for k in range(break_count):
        start = breaks[k]
        end = breaks[(k + 1) % break_count]

        if end > start:
            segment = pts[start : end + 1]
        else:
            segment = pts[start:] + pts[: end + 1]

        # Preserve hard corner endpoints; spline only the span between them.
        if len(segment) < 2:
            continue

        edge = make_open_spline_edge(name, segment)

        if hasattr(edge, "Edges"):
            edges.extend(edge.Edges)
        else:
            edges.append(edge)

    wire = close_edges_if_needed(edges)
    return wire


def make_profile_wire(name, points):
    if STEP_CURVE_MODE == "polyline":
        return make_polyline_wire(points)

    if STEP_CURVE_MODE not in ("spline", "spline-approx"):
        raise RuntimeError(f"Unsupported STEP_CURVE_MODE: {{STEP_CURVE_MODE}}")

    return make_segmented_spline_wire(name, points)


def make_extruded_profile(name, points, height_mm):
    try:
        wire = make_profile_wire(name, points)
    except Exception as e:
        raise RuntimeError(f"{{name}}: failed to create profile wire") from e

    if not wire.isValid():
        raise RuntimeError(f"{{name}}: generated wire is invalid; likely self-intersection, bad spline fit, or duplicate/tiny edges")

    if not wire.isClosed():
        raise RuntimeError(f"{{name}}: generated wire is not closed")

    try:
        face = Part.Face(wire)
    except Exception as e:
        raise RuntimeError(f"{{name}}: failed to create face; likely non-simple closed profile") from e

    if not face.isValid():
        raise RuntimeError(f"{{name}}: generated face is invalid; likely self-intersecting or non-simple profile")

    solid = face.extrude(FreeCAD.Vector(0, 0, height_mm))

    if not solid.isValid():
        raise RuntimeError(f"{{name}}: generated solid is invalid")

    return solid


doc = FreeCAD.newDocument("spring_cam_from_csv")

outer_solid = make_extruded_profile("outer", read_points(OUTER_CSV), SPRING_HEIGHT_MM)
inner_solid = make_extruded_profile("inner cutter", read_points(INNER_CSV), SPRING_HEIGHT_MM + 0.2)
cam_solid = make_extruded_profile("cam", read_points(CAM_CSV), CAM_HEIGHT_MM)

# Slightly oversize and shift the cutter in Z so the boolean cleanly cuts through.
inner_solid.translate(FreeCAD.Vector(0, 0, -0.1))

spring_solid = outer_solid.cut(inner_solid)
spring_solid = spring_solid.removeSplitter()

if not spring_solid.isValid():
    raise RuntimeError("spring: solid is invalid after boolean cut")

if not cam_solid.isValid():
    raise RuntimeError("cam: solid is invalid")

spring_obj = doc.addObject("Part::Feature", "Spring")
spring_obj.Shape = spring_solid

cam_obj = doc.addObject("Part::Feature", "Cam")
cam_obj.Shape = cam_solid

doc.recompute()

Import.export([spring_obj, cam_obj], str(OUTPUT_STEP))

if SPRING_OUTPUT_STEP:
    Import.export([spring_obj], SPRING_OUTPUT_STEP)

if CAM_OUTPUT_STEP:
    Import.export([cam_obj], CAM_OUTPUT_STEP)

print(f"Wrote {{OUTPUT_STEP}}")
if SPRING_OUTPUT_STEP:
    print(f"Wrote {{SPRING_OUTPUT_STEP}}")
if CAM_OUTPUT_STEP:
    print(f"Wrote {{CAM_OUTPUT_STEP}}")
'''
    path.write_text(freecad_script)


def export_step_with_freecad(
    freecadcmd: str,
    outer_csv: Path,
    inner_csv: Path,
    cam_csv: Path,
    output_step: Path,
    spring_output_step: Path | None,
    cam_output_step: Path | None,
    spring_height_mm: float,
    cam_height_mm: float,
    keep_freecad_script: bool,
    step_curve_mode: str,
    spline_approx_tolerance_mm: float,
    spline_break_angle_deg: float,
) -> None:
    output_step.parent.mkdir(parents=True, exist_ok=True)

    if keep_freecad_script:
        script_path = output_step.with_suffix('.freecad_export.py')
        write_freecad_export_script(
            script_path,
            outer_csv,
            inner_csv,
            cam_csv,
            output_step,
            spring_output_step,
            cam_output_step,
            spring_height_mm,
            cam_height_mm,
            step_curve_mode,
            spline_approx_tolerance_mm,
            spline_break_angle_deg,
        )
        subprocess.run([freecadcmd, str(script_path)], check=True)
        return

    with tempfile.TemporaryDirectory() as tmp:
        script_path = Path(tmp) / 'export_step.py'
        write_freecad_export_script(
            script_path,
            outer_csv,
            inner_csv,
            cam_csv,
            output_step,
            spring_output_step,
            cam_output_step,
            spring_height_mm,
            cam_height_mm,
            step_curve_mode,
            spline_approx_tolerance_mm,
            spline_break_angle_deg,
        )
        subprocess.run([freecadcmd, str(script_path)], check=True)


def resolve_cam_points(
    input_dir: Path,
    cam_choice: str,
    repair_cam: str,
) -> tuple[list[Point], str]:
    original_path = input_dir / 'cam_profile_m.csv'
    repaired_input_path = input_dir / 'cam_profile_repaired_m.csv'

    if cam_choice == 'repaired':
        cam_points = read_points_csv(repaired_input_path)
        return cam_points, f'existing repaired: {repaired_input_path.name}'

    if cam_choice == 'auto' and repaired_input_path.exists():
        cam_points = read_points_csv(repaired_input_path)
        return cam_points, f'existing repaired: {repaired_input_path.name}'

    cam_points = read_points_csv(original_path)

    intersections = count_self_intersections(cam_points, max_count=1)
    should_repair = repair_cam == 'on' or (repair_cam == 'auto' and intersections > 0)

    if should_repair:
        repaired_points, repair_count = repair_local_self_intersections(cam_points)
        remaining = count_self_intersections(repaired_points, max_count=1)

        if remaining:
            raise RuntimeError('cam repair did not remove all self-intersections')

        return repaired_points, f'auto-repaired from {original_path.name}; removed {repair_count} loops'

    return cam_points, f'original unrepaired: {original_path.name}'


def main() -> int:
    parser = argparse.ArgumentParser(description='Generate simplified OpenSCAD and optional STEP from spring-design-tool CSV files.')
    parser.add_argument('--input-dir', type=Path, default=Path.cwd(), help='Directory containing the CSV files.')
    parser.add_argument('--out-dir', type=Path, default=Path.cwd() / 'generated_spring_cam', help='Output directory.')
    parser.add_argument('--outer', default='outer_m.csv', help='Outer profile CSV filename.')
    parser.add_argument('--inner', default='inner_m.csv', help='Inner profile CSV filename.')
    parser.add_argument('--cam', choices=('auto', 'original', 'repaired'), default='auto', help='Cam profile source.')
    parser.add_argument('--repair-cam', choices=('auto', 'on', 'off'), default='auto', help='Repair cam self-intersections. Default: auto.')
    parser.add_argument('--simplify-tolerance', type=float, default=0.01, help='Polyline simplification tolerance in mm before STEP/SCAD generation. Use 0 to disable. Default: 0.01.')
    parser.add_argument('--step-curve-mode', choices=('polyline', 'spline', 'spline-approx'), default='spline', help='Curve type used only for STEP export. Default: spline.')
    parser.add_argument('--spline-approx-tolerance', type=float, default=0.02, help='Approximation tolerance in mm for --step-curve-mode spline-approx. Default: 0.02.')
    parser.add_argument('--spline-break-angle', type=float, default=25.0, help='Angle in degrees used to split spline spans at sharp corners. Default: 25.')
    parser.add_argument('--spring-height', type=float, default=3.0, help='Spring extrusion height in mm.')
    parser.add_argument('--cam-height', type=float, default=3.0, help='Cam extrusion height in mm.')
    parser.add_argument('--openscad', default=None, help='Path to openscad executable.')
    parser.add_argument('--freecadcmd', default=None, help='Path to FreeCADCmd/freecadcmd executable.')
    parser.add_argument('--no-open', action='store_true', help='Do not open the generated .scad file in OpenSCAD.')
    parser.add_argument('--step', action='store_true', default=True, help='Export STEP using FreeCADCmd if available. Default: enabled.')
    parser.add_argument('--no-step', action='store_false', dest='step', help='Disable STEP export.')
    parser.add_argument('--separate-step', action='store_true', help='Also export spring.step and cam.step separately.')
    parser.add_argument('--mesh', choices=('none', 'stl', '3mf'), default='none', help='Optional OpenSCAD mesh export format.')
    parser.add_argument('--mesh-mode', choices=('assembly', 'spring', 'cam'), default='assembly', help='Which body/bodies to export as mesh.')
    parser.add_argument('--keep-freecad-script', action='store_true', help='Keep the generated FreeCAD export script next to the STEP file.')

    args = parser.parse_args()

    if args.simplify_tolerance < 0.0:
        raise ValueError('--simplify-tolerance must be >= 0')

    if args.spline_approx_tolerance < 0.0:
        raise ValueError('--spline-approx-tolerance must be >= 0')

    if args.spline_break_angle <= 0.0 or args.spline_break_angle >= 180.0:
        raise ValueError('--spline-break-angle must be greater than 0 and less than 180')

    input_dir = args.input_dir.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    outer_csv = input_dir / args.outer
    inner_csv = input_dir / args.inner

    for path in (outer_csv, inner_csv):
        if not path.exists():
            raise FileNotFoundError(path)

    outer_points_original = read_points_csv(outer_csv)
    inner_points_original = read_points_csv(inner_csv)
    cam_points_original, cam_status = resolve_cam_points(
        input_dir=input_dir,
        cam_choice=args.cam,
        repair_cam=args.repair_cam,
    )

    outer_points = simplify_and_report('outer', outer_points_original, args.simplify_tolerance)
    inner_points = simplify_and_report('inner', inner_points_original, args.simplify_tolerance)
    cam_points = simplify_and_report('cam', cam_points_original, args.simplify_tolerance)

    simplified_outer_csv = out_dir / 'outer_simplified_m.csv'
    simplified_inner_csv = out_dir / 'inner_simplified_m.csv'
    simplified_cam_csv = out_dir / 'cam_simplified_m.csv'

    write_points_csv_mm_as_m(simplified_outer_csv, outer_points)
    write_points_csv_mm_as_m(simplified_inner_csv, inner_points)
    write_points_csv_mm_as_m(simplified_cam_csv, cam_points)

    scad_path = out_dir / 'spring_with_cam.scad'
    write_scad(
        scad_path,
        outer_points,
        inner_points,
        cam_points,
        args.spring_height,
        args.cam_height,
    )

    print(f'Wrote OpenSCAD file: {scad_path}')
    print(f'Wrote simplified CSV: {simplified_outer_csv}')
    print(f'Wrote simplified CSV: {simplified_inner_csv}')
    print(f'Wrote simplified CSV: {simplified_cam_csv}')
    print(f'Cam status: {cam_status}')
    print(f'STEP curve mode: {args.step_curve_mode}')

    openscad = find_executable(args.openscad, ('openscad', 'OpenSCAD'))
    freecadcmd = find_executable(args.freecadcmd, ('FreeCADCmd', 'freecadcmd', 'FreeCADCmd.exe', 'freecadcmd.exe'))

    if args.mesh != 'none':
        if not openscad:
            print('OpenSCAD executable not found; skipping mesh export.', file=sys.stderr)
        else:
            mesh_path = out_dir / f'spring_cam_{args.mesh_mode}.{args.mesh}'
            export_mesh_with_openscad(openscad, scad_path, mesh_path, args.mesh_mode)
            print(f'Wrote mesh: {mesh_path}')

    if args.step:
        if not freecadcmd:
            print('FreeCADCmd/freecadcmd not found; skipping STEP export.', file=sys.stderr)
            print('Install FreeCAD or pass --freecadcmd /path/to/FreeCADCmd.', file=sys.stderr)
        else:
            step_path = out_dir / 'spring_cam.step'
            spring_step_path = out_dir / 'spring.step' if args.separate_step else None
            cam_step_path = out_dir / 'cam.step' if args.separate_step else None

            export_step_with_freecad(
                freecadcmd=freecadcmd,
                outer_csv=simplified_outer_csv,
                inner_csv=simplified_inner_csv,
                cam_csv=simplified_cam_csv,
                output_step=step_path,
                spring_output_step=spring_step_path,
                cam_output_step=cam_step_path,
                spring_height_mm=args.spring_height,
                cam_height_mm=args.cam_height,
                keep_freecad_script=args.keep_freecad_script,
                step_curve_mode=args.step_curve_mode,
                spline_approx_tolerance_mm=args.spline_approx_tolerance,
                spline_break_angle_deg=args.spline_break_angle,
            )

    if not args.no_open:
        if not openscad:
            print('OpenSCAD executable not found; not opening GUI.', file=sys.stderr)
            print('Install OpenSCAD or pass --openscad /path/to/openscad.', file=sys.stderr)
        else:
            open_scad_gui(openscad, scad_path)
            print('Opened OpenSCAD GUI.')

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
