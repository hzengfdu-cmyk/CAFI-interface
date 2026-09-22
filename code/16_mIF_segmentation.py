"""Example Cellpose parameters for mIF segmentation."""

from pathlib import Path
import numpy as np
import tifffile
from cellpose import models
from skimage.measure import regionprops

input_image = Path("data/TMA_multichannel.tif")
output_dir = Path("results/cellpose")
output_dir.mkdir(parents=True, exist_ok=True)
tile_size = 1024
padding = 128

# The input array is CYX; update channel indices for the local image
channel = {"DAPI": 0, "EPCAM": 1, "CD8": 2, "PDGFRB": 3}
image = tifffile.imread(input_image).astype(np.float32)

def scale_channel(x, upper=99):
    x = np.maximum(x - np.percentile(x, 1), 0)
    low, high = np.percentile(x, [1, upper])
    return np.clip((x - low) / (high - low), 0, 1)

dapi = scale_channel(image[channel["DAPI"]], upper=99.8)
epcam = scale_channel(image[channel["EPCAM"]])
cd8 = scale_channel(image[channel["CD8"]])
pdgfrb = scale_channel(image[channel["PDGFRB"]])
composite = (epcam + cd8 + pdgfrb) / 3
segmentation_input = np.stack([composite, dapi], axis=-1)

model = models.CellposeModel(gpu=True, pretrained_model="cpsam")
# The full analysis used 1,024 × 1,024 core tiles with 128-pixel context on each side.
# This script shows a single-tile call; tiled processing retains objects centred in the core tile.
masks, _, _ = model.eval(
    segmentation_input,
    channel_axis=-1,
    normalize=True,
    diameter=30,
    flow_threshold=0.4,
    cellprob_threshold=0,
    min_size=80,
    tile_overlap=0.1,
)
tifffile.imwrite(output_dir / "cell_masks.tif", masks.astype(np.uint32))

# Extract centroids, areas and mean intensities from the original channels
rows = []
for cell in regionprops(masks):
    y, x = cell.centroid
    rr, cc = cell.coords.T
    row = {"cell_id": cell.label, "centroid_x": x, "centroid_y": y, "area_px": cell.area}
    for marker, idx in channel.items():
        row[marker] = float(image[idx, rr, cc].mean())
    rows.append(row)

import pandas as pd
pd.DataFrame(rows).to_csv(output_dir / "cell_table.tsv", sep="\t", index=False)
