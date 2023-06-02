# BloodDose
BloodDose: a revised version of HEDOS. Blood flow is stochastically modeled and superimposed with a - potentially time varying - dose rate to calculate blood dose.

#
For blood dose calculations, the configuration parameters (patient, treatment and simulation parameters) are set up in **BloodDose.py**.
It then calls one of two possible workflows: 
- **BloodDoseFromFields.py** calculates the blood dose given the patient's dose distribution and the segmentations of the organs which are contributing to the overall blood dose. For the example here we have used a mesh reference phantom (ICRP Publication 145) and create an artificial sample dose. For patient-specific blood dose calculations, this should be replaced with patient data.
- **BloodDoseFromDVHs.py** does the same thing, but uses DVHs of each of the organs contributing to the overall blood dose.

#
The calculation of blood dose follows these steps in succession:
- **FlowModel.py**: Set up a graph that reflects the connectivity and magnitude of blood flow between blood compartments. Convert this into a matrix of transition probabilities.  
- **TemporalDistribution.py**: Simulate the blood flow over time. Blood particles flow through the model by a stochastic jumping process based on survival analysis.
- **CompartmentDose.py**: Accumulate dose in the blood particles.


