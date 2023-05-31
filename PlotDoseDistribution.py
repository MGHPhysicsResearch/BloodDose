import numpy as np
import matplotlib.pyplot as plt


def plot_dose_distribution(blood_dose_total, dose_contributions, mean_blood_dose=None):
    x_max1 = 2 * np.percentile(blood_dose_total.dose, 90)
    x_max2 = max([2 * np.percentile(dose, 90) for dose in dose_contributions.values()])
    bins1 = np.linspace(0, x_max1, 100)
    bins2 = np.linspace(0, x_max2, 100)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.hist(blood_dose_total.dose, bins=bins1)
    ax1.axvline(np.mean(blood_dose_total.dose), ymax=1, c='k', linestyle='-',
                label='Simulated mean dose - {:.3f} Gy'.format(np.mean(blood_dose_total.dose)))
    if mean_blood_dose is not None:
        ax1.axvline(mean_blood_dose, ymax=1, c='red', linestyle='--',
                    label='Expected mean dose - {:.3f} Gy'.format(mean_blood_dose))
    for organ, dose in dose_contributions.items():
        ax2.hist(dose, bins=bins2[1:], histtype='step', linewidth=1.5, alpha=0.7, label=organ)

    ax1.set_xlabel('Dose (Gy)')
    ax2.set_xlabel('Dose (Gy)')
    ax1.set_title('Blood dose histogram')
    ax2.set_title('Blood dose contributions')
    ax1.legend()
    ax2.legend()

    plt.show()

