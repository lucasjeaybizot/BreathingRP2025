# BreathingRP2025
This repository contains the code for Jeay-Bizot et al. (2025) Matters Arising following Park et al. (2020) article.



### Repository Structure

````
.
├── binning\_corrected.m
├── binning\_original.m
├── generate\_simulateddata.m
├── make\_figure1.m
├── make\_figureS1.m ─ make\_figureS5.m
├── preprocess\_rawdata.m
├── src/
├── figures/               # Output figures
├── results/
└── data/                  # needs to be manually added by running make_datadirs.m

````

---

## How to Run the Code

### 1. Preprocessing

Start by preprocessing the raw EEG and respiration signals:

```matlab
preprocess_rawdata.m
````

This step expects the `data/new/raw/` folder to be present.

### 2. Binning Trials by Respiration Phase

Run both of the following scripts to bin the data (saved to results\):

```matlab
binning_corrected.m
```

and

```matlab
binning_original.m
```

This step expects all the `data/.../processed/` folder to be present.

### 3. Generate Figures

Run the figure scripts in the following order:

```matlab
make_figure1.m
make_figureS1.m
make_figureS2.m
make_figureS3.m
make_figureS4.m
make_figureS5.m
```

Each script generates one of the figures presented in the manuscript.

---

## Data Availability

The dataset used in this project consists of:

* **`data/new/`**: Raw and processed data collected for this study (**publicly available**)
* **`data/simulated/`**: Synthetic data used for validation (**run generate_simulateddata.m**)
* **`data/original/`**: Raw data from Park et al. (2020) — **not included**

### Download Public Data (\~10 GB)

You can download the `data/new/raw` from our OSF project:

[https://osf.io/3cvyk/?view_only=e02aa800fc5d40a082ebd2c862f768c0](https://osf.io/3cvyk/?view_only=e02aa800fc5d40a082ebd2c862f768c0)](https://osf.io/3cvyk/)

### Accessing Restricted Data

To obtain the original dataset from Park et al. (2020), please contact the authors of that study directly.

### Generating simulated data

Run the generate_simulateddata.m script to populate data\simulated\processed.

---

## Requirements

* MATLAB R2020b or later
* FieldTrip toolbox (latest version from [https://www.fieldtriptoolbox.org](https://www.fieldtriptoolbox.org))
* circ_stat toolbox (latest version from [https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics](https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics))

---

## Notes

* For questions or feedback, contact: [jeaybizot@chapman.edu](mailto:jeaybizot@chapman.edu)

---

