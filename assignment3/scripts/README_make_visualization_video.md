# make_visualization_video.py

Create one side-by-side video from screenshot subfolders under `assignment3/visualization/<run-folder>`.

## Configure in script

Edit the user variables at the top of `assignment3/scripts/make_visualization_video.py`:

- `RUN_FOLDER`
- `FIELDS`
- `FPS`
- `OUTPUT`
- `MAX_FRAMES`

Supported fields (aliases accepted):

- `Cell_Field`
- `Oxygen`
- `VEGF2` (or `VEGF 2`)

## Run

```powershell
python assignment3/scripts/make_visualization_video.py
```

Default output (if `OUTPUT = None`):

- `assignment3/visualization/<run-folder>/combined_<fields>.mp4`
- If MP4 backend is unavailable, the script automatically writes a `.gif` instead.

The script aligns panels by common timestep (intersection across selected fields) and places panels side-by-side in the same order as `FIELDS`.
