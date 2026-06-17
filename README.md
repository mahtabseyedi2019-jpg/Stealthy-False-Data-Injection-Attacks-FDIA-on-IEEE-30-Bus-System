# Stealthy-False-Data-Injection-Attacks-FDIA-on-IEEE-30-Bus-System

Abstract
Attacks (FDIA) engineered for the standard IEEE 30-bus power system. Unlike conventional datasets that focus on easily detectable anomalies, this dataset emphasizes highly stealthy micro-injections mathematically designed to bypass standard Bad Data Detection (BDD) mechanisms, specifically the Largest Normalized Residual (LNR) and Chi-Squared tests.

The baseline normal operations are generated using Direct Current Power Flow (DCPF) equations integrated with real-world, highly variable load profiles sourced from the New York Independent System Operator (NYISO). To evaluate the robustness of modern machine learning and deep learning architectures, the attacks are categorized into five progressive scenarios. These scenarios range from High-Intensity attacks (Scenario 1) to an Advanced Combined scenario (Scenario 5), where the attack magnitudes are strictly restricted to a highly subtle range of [0.005, 0.05] radians, resulting in physically imperceptible angular deviations.

To maximize usability, the dataset is provided in two distinct formats:

CSV Format: 2D tabular data containing branch active power flows and multi-label attack statuses, ideal for statistical analysis and classical machine learning.
NPY Format: 3D time-series arrays (Samples, Timesteps, Features) pre-processed with a sliding window length of 16, feature scaling, and class balancing (1:2 ratio), making it directly deployable for deep learning architectures (e.g., LSTMs, CNNs) targeting the dual challenges of attack Detection and Localization.
Instructions:
The dataset is organized into 5 progressive scenarios. Each scenario is provided in two formats: 1. CSV Format: 2D tabular data for classical ML and statistical analysis. 2. NPY Format: 3D time-series arrays pre-formatted for deep learning models (e.g., LSTMs, CNNs).

Preprocessing Details (NPY files):

Temporal Windowing: A sliding window of 16 timesteps is applied.
Feature Scaling: Standard scaling is applied (fitted strictly on the training set to prevent data leakage).
Class Balancing: The normal (benign) class is undersampled to achieve a 1:2 ratio (Attacked vs. Normal) to prevent class dominance during training.
Suggested AI Applications:

Binary Classification for rapid anomaly detection.
Multi-label Classification for precise attack localization across buses.
Benchmarking model robustness against decreasing attack intensities (from Scenario 1 to 5).
