import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def vol_to_gridpoints(vol, affine):
    dims = np.array(vol.shape)
    center_voxel = (dims + 1) / 2
    center = np.dot(affine, np.append(center_voxel, 1))[:3]
    extents = np.dot(affine[:3, :3], dims)
    signs = np.sign(np.diag(affine[:3, :3]))

    # Create the gridpoints at the center of the voxels:
    x = np.linspace(center[0] - 0.5 * extents[0], center[0] + 0.5 * extents[0], dims[0]) * signs[0]
    y = np.linspace(center[1] - 0.5 * extents[1], center[1] + 0.5 * extents[1], dims[1]) * signs[1]
    z = np.linspace(center[2] - 0.5 * extents[2], center[2] + 0.5 * extents[2], dims[2]) * signs[2]
    return x, y, z


def plot_volumes(volume_ref, volume, plot_slice=None, cmap_ref='jet', cmap='Greys_r'):
    if plot_slice is None:
        plot_slice = np.argmax(np.sum(volume, axis=(0, 1)))
    plt.imshow(volume_ref[:, :, plot_slice], cmap=cmap_ref)
    c_bar = plt.colorbar()
    c_bar.set_label('Treatment dose at slice {} (Gy)'.format(plot_slice))
    plt.imshow(volume[:, :, plot_slice], cmap=cmap, alpha=0.5)
    plt.show()


class Patient:
    """
    Here we load the dose and segmentations simply as numpy files.
    Hence this requires that you have resampled everything to the same grid (e.g. that of the patient scan first).
    Furthermore, as numpy arrays are just arrays (no coordinates), we also load the affine transform that positions it
    in space and convert the grid to x, y, z coordinates (the gridpoints).
    """
    def __init__(self):
        self.gridpoints = None
        self.dose = None
        self.tumor_volume_fraction = None
        self.seg_organs = {}

    def read_from_numpy(self, read_dir, organ_names, plot=True):
        """
        Read in segmentations which are assumed to be bundled in a .npz files.
        Read in dose (which has been created artificially, just as an example).
        Read in an affine transform which defines the coordinates of the voxels of the numpy arrays.

        This could/should be replaced by your own function, potentially reading in DICOM files of patients directly.
        """
        segs_loaded = np.load(os.path.join(read_dir, 'compressed_segs.npz'))
        self.dose = np.load(os.path.join(read_dir, 'dose.npy'))
        affine = np.load(os.path.join(read_dir, 'affine.npy'))
        self.gridpoints = vol_to_gridpoints(self.dose, affine)
        for organ_name in organ_names:
            self.seg_organs[organ_name] = segs_loaded[organ_name]
            if plot:
                plot_volumes(self.dose, self.seg_organs[organ_name])

    def get_tumor_volume_fraction(self, tumor_bearing_organ, tumor):
        """
        This function gets the volume fraction of the tumor with respect to the organ in which it resides.
        """
        self.tumor_volume_fraction = np.sum(self.seg_organs[tumor]) / np.sum(self.seg_organs[tumor_bearing_organ])
        print('Tumor volume fraction = {:.4f}'.format(self.tumor_volume_fraction))

    def get_mean_organ_dose(self, organ_name):
        """
        This function calculates the mean organ dose.
        """
        idx = np.where(self.seg_organs[organ_name] == 1)
        mod = np.mean(self.dose[idx[0], idx[1], idx[2]])
        return mod

    def write_dvh(self, save_dir, organ_names):
        """
        Summarize fields into DVHs of each organ separately. Write out in csv-files.
        This representation can also be used for blood dose calculation.
        """
        bins = np.arange(0, np.ceil(np.max(self.dose)) + 0.1, 0.1)
        for organ_name in organ_names:
            idx = np.where(self.seg_organs[organ_name] == 1)
            organ_dose = self.dose[idx[0], idx[1], idx[2]]
            values, bins = np.histogram(organ_dose, bins=bins)
            coverage = np.append(np.cumsum(values[::-1])[::-1] / organ_dose.size * 100, [0])
            dvh = np.stack([bins, coverage], axis=1)
            # save
            pd.DataFrame(dvh).to_csv(os.path.join(save_dir, organ_name + '_DVH.csv'), header=['dose_bins', 'coverage'])
            plt.plot(bins, coverage, label=organ_name)
        plt.legend()
        plt.xlabel('Dose (Gy)')
        plt.ylabel('Coverage (%)')
        plt.show()

