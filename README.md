# Compensating Internal Model Control for Power Inverters

This repository contains simulation and research material developed for a thesis on a new approach to inverter current control based on compensating internal models.

The work focuses on current regulation for inverter systems in the presence of periodic references and disturbances. The main objective is to study a control framework that is better suited for implementation and adaptation than standard internal-model-based designs, while still preserving the essential regulation properties required for accurate tracking and disturbance rejection.

The repository includes MATLAB scripts, Simulink models, thesis documents, and supporting third-party code used in the development and validation of the proposed approach.

## Overview

The main topic of this work is a compensating internal model approach for regulator design, with application to inverter current control.

In inverter control problems, periodic signals and harmonics are fundamental. Classical approaches such as repetitive control and resonant control are widely used for this purpose. This thesis investigates an alternative formulation based on compensating internal models, with particular attention to:

- adaptation of internal model parameters
- comparison with the H-infinity repetitive control method
- application to inverter systems with LCL filters
- practical simulation-based validation in MATLAB and Simulink

The repository is therefore centered on a new current-control approach that is mainly based on the theory developed in the papers listed below, together with the simulation models used in the thesis.

## Main idea

The central idea of the thesis is to use compensating internal models for inverter current control.

Rather than treating the internal model only as a fixed embedded block inside a classical regulator, this approach studies a more structured compensation framework that can be better adapted to implementation requirements and, when needed, extended with parameter adaptation.

This is relevant for inverter applications because:

- current references are often sinusoidal or periodic
- harmonic compensation is essential
- robustness to model uncertainty is important
- adaptation can improve performance when frequencies or operating conditions vary

The work is therefore not limited to standard resonant or repetitive designs, but investigates a more general regulator framework and its application to inverter control problems.

## Contents of the repository

## LC filter

The `LC filter` folder contains models and initialization scripts illustrating the basic idea of state-feedback control with compensation of the internal model and parameter adaptation.

This part is not an inverter-bridge model. It is a simpler control setup built around an LC filter, disturbance input, and control signal, and is used to demonstrate the core concept of the method before moving to more application-oriented cases.

## Repetitive Inverter Control

The `Repetitive Inverter Control` folder contains models connected to repetitive-control-type structures and comparisons.

This part is useful for studying the relation between the proposed compensating internal model approach and more classical repetitive control methods for periodic signal compensation.

The folder includes:

- CIM-based repetitive control models
- adaptive variants
- an H-infinity comparison model

## Resonant Inverter Control

The `Resonant Inverter Control` folder contains proportional-resonant-type control models combined with the compensating internal model framework.

Included cases:

- adaptive, single-output control
- non-adaptive, single-output control
- non-adaptive, two-output control

These models are relevant because resonant control is a standard tool in inverter current regulation, especially in applications involving sinusoidal references and harmonic rejection.

## Thesis scope

The thesis studies inverter current control as a regulator problem and develops a framework based on compensating internal models.

The main goals are:

- to formulate inverter current regulation problems in a regulator-theoretic setting
- to apply compensating internal model structures to periodic tracking and disturbance rejection
- to study implementation-oriented controller forms
- to investigate adaptive extensions
- to compare the proposed method with repetitive alternatives
- to validate the approach through MATLAB and Simulink simulations

The two PDF files in the root folder contain the thesis manuscript and the executive summary and provide the full theoretical background, methodology, derivations, and results.

## Main references

This work is mainly based on the following publications:

- P. Colaneri, G. P. Incremona, L. Marconi, and L. Mirkin. On the implementation and adaptation of a class of internal models. In 25th International Symposium on Mathematical Theory of Networks and Systems MTNS 2022, pages 968–971, Bayreuth, Germany, 2022.

- P. Colaneri, G. P. Incremona, and L. Mirkin. On internal model compensation for general regulator problems. In 62nd IEEE Conference on the Decision and Control (CDC), pages 2657–2662, 2023. doi: 10.1109/CDC49753.2023.10383453.

- P. Colaneri, G. P. Incremona, and L. Mirkin. On compensating internal models in regulator problems. IEEE Transactions on Automatic Control, pages 1–13, 2025. doi: 10.1109/TAC.2025.3606835.

- G. P. Incremona, L. Mirkin, and P. Colaneri. Integral sliding-mode control with internal model: A separation. IEEE Control Systems Letters, 6:446–451, 2021. doi: 10.1109/LCSYS.2021.3079187.

- G. P. Incremona, P. Colaneri, and L. Mirkin. Parameter adaptation for general regulator problems with compensation of internal model. In 63rd IEEE Conference on the Decision and Control (CDC), pages 2934–2939, 2024. doi: 10.1109/CDC56724.2024.10886230.

- L. Mirkin. On dead-time compensation in repetitive control. IEEE Control Systems Letters, 4(4):791–796, 2020. doi: 10.1109/LCSYS.2020.2992712.

## Third-party code

The repository includes third-party material in:

`third_party/generalized_sylvester/`

This folder contains:

- `generalized_sylvester.m`
- `LICENSE.txt`

This code is included as an external dependency and should be treated according to its own license terms. It is not part of the original thesis contribution itself, but is used as supporting material within the implementation workflow.

## Software environment

The repository is intended for use with MATLAB and Simulink.

Main file types:

- `.m` files for initialization and parameter setup
- `.slx` files for simulation models
- `.pdf` files for the thesis and executive summary

To reproduce results, open the relevant initialization script first and then run the corresponding Simulink model.

## Notes

This repository is primarily a research and thesis repository. Its purpose is to document and reproduce the simulation framework used to study compensating internal model control for inverter applications.

It combines:

- theoretical material from the thesis
- simulation models used for validation
- comparisons with alternative controller structures
- third-party support code required for some computations
