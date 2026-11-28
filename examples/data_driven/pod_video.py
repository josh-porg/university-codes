"""POD of video frames: low-order reconstructions of an ignition / wavy-surface video
(``PODwavy``, ``Proper_orthogonal_decomposition_Grey_Ignition_video``,
``Proper_orthogonal_decomposition_RGB_Ignition``, ``imageVisualizationTest``,
``video2array``, ``SaveVideoToArray``).

Input is a ``.mat`` holding a frame array (``videoFrames`` H x W x 3 x T or
``greyVideoFrames``) or a video file read with ``imageio`` (optional
dependency)::

    python pod_video.py wavy_video.mat --step 10 --modes 7 --out wavy_pod.gif
    python pod_video.py ignition.mp4 --grey --modes 5

The reconstruction adds modes one at a time to the mean as in the MATLAB
loops; ``--out`` writes the side-by-side comparison as an animation.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from unicodes.decomposition import POD
from unicodes.io import load_mat


def video_to_array(path):
    """Frames of a video as ``(H, W, 3, T)`` (``video2array``)."""
    import imageio.v3 as iio

    return np.moveaxis(np.asarray(iio.imread(path)), 0, -1)


p = argparse.ArgumentParser()
p.add_argument("input")
p.add_argument("--variable", default=None)
p.add_argument("--step", type=int, default=10, help="spatial subsampling")
p.add_argument("--grey", action="store_true", help="use the first colour channel only")
p.add_argument("--modes", type=int, default=5)
p.add_argument("--frames", type=int, default=55)
p.add_argument("--out", default=None)
a = p.parse_args()

if a.input.endswith(".mat"):
    d = load_mat(a.input)
    key = a.variable or next(k for k in ("videoFrames", "greyVideoFrames") if k in d)
    frames = np.asarray(d[key])
else:
    frames = video_to_array(a.input)
frames = frames[:: a.step, :: a.step].astype(float)
if a.grey and frames.ndim == 4:
    frames = frames[:, :, 0, :]
pod = POD(frames)
print("energy fraction of the first modes:", np.round(pod.energy_fraction[: a.modes + 3], 4))
recon = np.clip(pod.reconstruct(a.modes), 0, 255)


def show(ax, img):
    ax.clear()
    ax.imshow(img.astype(np.uint8), cmap=None if img.ndim == 3 else "gray", vmin=0, vmax=255)
    ax.axis("off")


fig, (ax1, ax2) = plt.subplots(1, 2)


def frame(n):
    show(ax1, recon[..., n])
    ax1.set_title(f"mean + {a.modes} modes")
    show(ax2, frames[..., n])
    ax2.set_title("raw data")


frame(0)
if a.out:
    FuncAnimation(fig, frame, frames=min(a.frames, frames.shape[-1]), interval=50).save(a.out)
plt.show()
