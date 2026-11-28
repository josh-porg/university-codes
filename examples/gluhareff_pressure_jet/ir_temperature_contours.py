"""R. Bramlette (University of Kansas Ph.D. research), ``Bramlette_ku_0099D_14625_DATA_2IR`` (IMG2IR):
apparent-temperature contours from a photograph of a part glowing at incandescence.

The image is converted to grey levels (MATLAB ``rgb2gray`` weights) and the
brightness is mapped linearly onto the measured temperature range
(default 800-2000 deg F), then contoured every ``--step`` degrees. The
MATLAB scaled ``BW (Tmax - Tmin) / (Imax - Imin) + Tmin`` without
subtracting ``Imin``, which shifts every temperature up unless the darkest
pixel is black; ``(BW - Imin)`` is used here. The file dialog is the
``image`` argument::

    python ir_temperature_contours.py engine.png --t-max 2000 --t-min 800 --step 10
"""

import argparse
import time

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

p = argparse.ArgumentParser()
p.add_argument("image")
p.add_argument("--t-max", type=float, default=2000.0, help="maximum temperature, deg F")
p.add_argument("--t-min", type=float, default=800.0, help="minimum temperature, deg F")
p.add_argument("--step", type=float, default=10.0, help="contour step, deg F")
a = p.parse_args()

t0 = time.process_time()
img = np.asarray(Image.open(a.image).convert("RGB"), float)
rows, cols = img.shape[:2]
print(f"{a.image} was read into memory ({rows} pixels tall, {cols} pixels wide)")
bw = np.round(img @ np.array([0.2989, 0.5870, 0.1140]))  # rgb2gray, uint8 levels
i_max, i_min = bw.max(), bw.min()
T_app = (bw - i_min) * (a.t_max - a.t_min) / (i_max - i_min) + a.t_min

fig, ax = plt.subplots(2, 1, figsize=(7, 8))
ax[0].imshow(img.astype(np.uint8))
ax[0].set_title(f"Original image: {a.image}")
ax[0].axis("off")
print("Generating contours of apparent temperature from the image file...")
levels = np.arange(a.t_min, a.t_max + a.step / 2, a.step)
cf = ax[1].contourf(T_app, levels, cmap="inferno")
if len(levels) <= 10:
    ax[1].contour(T_app, levels, colors="k", linewidths=0.5)
ax[1].set(xlim=(0, cols), ylim=(rows, 0), aspect="equal", xticks=[], yticks=[],
          title=f"Apparent temperature contours, $\\Delta T$ = {a.step:.0f} °F")
fig.colorbar(cf, ax=ax[1], orientation="horizontal", label="Apparent temperature (°F)")
print(f"Contour image processing complete in {time.process_time() - t0:.1f} s")
plt.show()
