import re
from pathlib import Path

import imageio.v2 as imageio
import numpy as np
from PIL import Image, ImageDraw, ImageFont


FIELD_ALIASES = {
    "cell_field": "Cell_Field",
    "cell field": "Cell_Field",
    "cell": "Cell_Field",
    "oxygen": "Oxygen",
    "vegf2": "VEGF2",
    "vegf 2": "VEGF2",
    "vegf": "VEGF2",
}

FIELD_SUBFOLDER_HINTS = {
    "Cell_Field": ["cell_field"],
    "Oxygen": ["oxygen"],
    "VEGF2": ["vegf2", "vegf_2", "vegf 2"],
}

# User configuration (edit these values and run the script)
RUN_FOLDER = "tumor_growth"
FIELDS = ["Cell_Field", "Oxygen", "VEGF2"]
FPS = 6
OUTPUT = "../visualization/tumor_growth/tumor_growth.mp4"
MAX_FRAMES = None  # Example: 200

# Panel label styling
PANEL_LABEL_HEIGHT = 44
PANEL_LABEL_FONT_SIZE = 24

# Per-panel labels (keys must use canonical field names)
PANEL_LABELS = {
    "Cell_Field": "Cell distribution",
    "Oxygen": "Oxygen concentration",
    "VEGF2": "VEGF concentration",
}

# Cell-field legend entries (user-editable).
# Each item requires: name + color (hex code).
CELL_FIELD_LEGEND = [
    {"name": "Medium", "color": "#000000"},
    {"name": "Normal", "color": "#0000ff"},
    {"name": "Hypoxic", "color": "#008000"},
    {"name": "Necrotic", "color": "#ff0000"},
    # {"name": "ActiveNeovascular", "color": "#DC1E1E"},
    {"name": "Vascular", "color": "#aa0000"},
    # {"name": "InactiveNeovascular", "color": "#B400B4"},
]


def canonical_field_name(raw_name: str) -> str:
    key = raw_name.strip().lower().replace("_", " ")
    canonical = FIELD_ALIASES.get(key)
    if canonical is None:
        raise ValueError(
            f"Unsupported field '{raw_name}'. Use one of: Cell_Field, Oxygen, VEGF2/VEGF 2"
        )
    return canonical


def resolve_visualization_root() -> Path:
    return (Path(__file__).resolve().parents[1] / "visualization").resolve()


def resolve_run_folder(root: Path, run_folder: str) -> Path:
    run_path = (root / run_folder).resolve()
    if not run_path.exists() or not run_path.is_dir():
        raise FileNotFoundError(f"Run folder not found: {run_path}")
    return run_path


def find_field_subfolder(run_path: Path, canonical_field: str) -> Path:
    hints = FIELD_SUBFOLDER_HINTS[canonical_field]
    candidates = []
    for child in run_path.iterdir():
        if not child.is_dir():
            continue
        lowered = child.name.lower()
        if any(hint in lowered for hint in hints):
            candidates.append(child)

    if not candidates:
        raise FileNotFoundError(
            f"No subfolder found for field '{canonical_field}' inside {run_path}"
        )

    candidates.sort(key=lambda p: p.name)
    return candidates[0]


def extract_timestep(image_path: Path) -> int:
    match = re.search(r"_(\d+)\.png$", image_path.name, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"Could not parse timestep from filename: {image_path.name}")
    return int(match.group(1))


def collect_frames(field_folder: Path) -> dict[int, Path]:
    pngs = sorted(field_folder.glob("*.png"))
    if not pngs:
        raise FileNotFoundError(f"No PNG frames in {field_folder}")

    by_time: dict[int, Path] = {}
    for path in pngs:
        timestep = extract_timestep(path)
        by_time[timestep] = path
    return by_time


def _hex_to_rgb(hex_code: str) -> tuple[int, int, int]:
    code = hex_code.strip().lstrip("#")
    if len(code) != 6 or not re.fullmatch(r"[0-9a-fA-F]{6}", code):
        raise ValueError(f"Invalid hex color '{hex_code}'. Expected format '#RRGGBB'.")
    return tuple(int(code[i : i + 2], 16) for i in (0, 2, 4))


def _ensure_rgb(frame: np.ndarray) -> np.ndarray:
    if frame.ndim == 2:
        return np.stack([frame, frame, frame], axis=-1)
    if frame.shape[2] == 4:
        return frame[:, :, :3]
    return frame


def _load_label_font() -> ImageFont.ImageFont:
    for font_name in ("arial.ttf", "DejaVuSans-Bold.ttf", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(font_name, PANEL_LABEL_FONT_SIZE)
        except Exception:
            continue
    return ImageFont.load_default()


def _add_panel_label(frame: np.ndarray, label: str) -> np.ndarray:
    label_height = int(max(28, PANEL_LABEL_HEIGHT))
    image = Image.new("RGB", (frame.shape[1], frame.shape[0] + label_height), (20, 20, 20))
    image.paste(Image.fromarray(frame), (0, label_height))

    draw = ImageDraw.Draw(image)
    font = _load_label_font()
    text = str(label)
    text_bbox = draw.textbbox((0, 0), text, font=font)
    text_w = text_bbox[2] - text_bbox[0]
    text_h = text_bbox[3] - text_bbox[1]

    text_x = max(0, (frame.shape[1] - text_w) // 2)
    text_y = max(0, (label_height - text_h) // 2)
    draw.text((text_x, text_y), text, fill=(255, 255, 255), font=font)
    return np.array(image)


def _add_cell_field_legend(frame: np.ndarray) -> np.ndarray:
    if not CELL_FIELD_LEGEND:
        return frame

    image = Image.fromarray(frame)
    draw = ImageDraw.Draw(image)
    font = ImageFont.load_default()

    swatch = 12
    row_h = 18
    legend_title = "Legend"
    content_w = max(
        [draw.textlength(legend_title, font=font)]
        + [draw.textlength(str(item["name"]), font=font) for item in CELL_FIELD_LEGEND]
    )
    box_w = int(28 + swatch + content_w)
    box_h = int(10 + row_h * (len(CELL_FIELD_LEGEND) + 1))

    x0 = max(0, frame.shape[1] - box_w - 8)
    y0 = 8
    x1 = min(frame.shape[1], x0 + box_w)
    y1 = min(frame.shape[0], y0 + box_h)

    draw.rectangle((x0, y0, x1, y1), fill=(15, 15, 15), outline=(220, 220, 220), width=1)
    draw.text((x0 + 8, y0 + 4), legend_title, fill=(255, 255, 255), font=font)

    for idx, item in enumerate(CELL_FIELD_LEGEND):
        color = _hex_to_rgb(str(item["color"]))
        name = str(item["name"])
        y = y0 + 4 + (idx + 1) * row_h
        draw.rectangle((x0 + 8, y, x0 + 8 + swatch, y + swatch), fill=color, outline=(240, 240, 240), width=1)
        draw.text((x0 + 8 + swatch + 6, y - 1), name, fill=(255, 255, 255), font=font)

    return np.array(image)


def pad_to_height(frame: np.ndarray, target_height: int) -> np.ndarray:
    if frame.shape[0] == target_height:
        return frame
    pad_rows = target_height - frame.shape[0]
    if pad_rows < 0:
        return frame[:target_height, :, :]
    return np.pad(frame, ((0, pad_rows), (0, 0), (0, 0)), mode="constant", constant_values=0)


def compose_frame_row(panel_specs: list[tuple[str, Path]]) -> np.ndarray:
    rgb_frames = []
    for field_name, frame_path in panel_specs:
        frame = _ensure_rgb(imageio.imread(str(frame_path)))
        if field_name == "Cell_Field":
            frame = _add_cell_field_legend(frame)
        frame = _add_panel_label(frame, PANEL_LABELS.get(field_name, field_name))
        rgb_frames.append(frame)

    max_height = max(frame.shape[0] for frame in rgb_frames)
    aligned = [pad_to_height(frame, max_height) for frame in rgb_frames]
    return np.concatenate(aligned, axis=1)


def default_output_path(run_path: Path, selected_fields: list[str]) -> Path:
    suffix = "_".join(field.lower() for field in selected_fields)
    return run_path / f"combined_{suffix}.mp4"


def open_writer_with_fallback(output_path: Path, fps: int):
    try:
        writer = imageio.get_writer(str(output_path), fps=fps, macro_block_size=None)
        return writer, output_path
    except Exception as exc:
        gif_path = output_path.with_suffix(".gif")
        print(
            "[make_visualization_video] MP4 backend unavailable; "
            f"falling back to GIF ({exc})."
        )
        writer = imageio.get_writer(str(gif_path), fps=fps)
        return writer, gif_path


def main() -> None:
    selected_fields = [canonical_field_name(name) for name in FIELDS]

    # Keep order but remove duplicates.
    deduped_fields = list(dict.fromkeys(selected_fields))

    viz_root = resolve_visualization_root()
    run_path = resolve_run_folder(viz_root, RUN_FOLDER)

    frame_maps = {}
    for field_name in deduped_fields:
        subfolder = find_field_subfolder(run_path, field_name)
        frame_maps[field_name] = collect_frames(subfolder)

    common_times = sorted(set.intersection(*(set(m.keys()) for m in frame_maps.values())))
    if not common_times:
        raise RuntimeError("No common timesteps across selected fields.")

    if MAX_FRAMES is not None:
        common_times = common_times[: max(0, int(MAX_FRAMES))]
    if not common_times:
        raise RuntimeError("No frames selected for output after applying MAX_FRAMES.")

    output_path = Path(OUTPUT).resolve() if OUTPUT else default_output_path(run_path, deduped_fields)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fps = max(1, int(FPS))
    writer, final_output_path = open_writer_with_fallback(output_path, fps)
    with writer:
        for timestep in common_times:
            panel_specs = [(field, frame_maps[field][timestep]) for field in deduped_fields]
            panel_frame = compose_frame_row(panel_specs)
            writer.append_data(panel_frame)

    print(f"[make_visualization_video] Saved {len(common_times)} frames to: {final_output_path}")
    print(f"[make_visualization_video] Fields: {', '.join(deduped_fields)}")


if __name__ == "__main__":
    main()

