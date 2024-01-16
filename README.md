# EngineeringReasoningCourse

## Background
As college readiness continues to decline, the proportion of students entering engineering programs with low mathematics proficiency is increasing. These students have lower retention rates than their calculus-ready peers. 
## Purpose
This work describes the design and implementation of an Introduction to Engineering Problem Solving course targeted to first-year engineering students who place into College Algebra with the goal of increasing student retention in engineering.
## Design/Method
81 students were enrolled in the course over two years. The course implements the CDIO educational framework and Paul-Elder Critical Thinking Theory to promote engineering problem solving and critical thinking skills through problem and project-based learning. Institutional data including cumulative GPA, mathematics course grades, and enrollment status and within the institution were collected for four semesters. The experimental group took the Critical Thinking Assessment Test (CAT) at the start and end of the intervention. All statistical analysis was completed in R using appropriate Bayesian statistical methods.
## Results
Students participating in the intervention course saw improved critical thinking skills compared to their baseline along with higher pass rates in their College Algebra course and improved progression through the mathematics curriculum compared to the control group. Cumulative GPA and retention within engineering and the university were also improved.
## Conclusion
The importance of supporting non-calculus ready engineering students will continue to increase in the coming years. Intervention coursework has the potential to significantly increase student success in engineering. Future work should explore the possibility of multi-semester interventions to support students until they reach calculus.

## Statistical Analysis
All statistical analysis was performed in R v4.3.1. Bayesian estimation using the BEST (Bayesian Estimation Supersedes the t Test) method developed by John K. Kruschke was completed using The BEST package (v0.5.4) to estimate the size of the effect the intervention course had on student parameters such as first-semester cumulative GPA (Kruschke, 2013). Comparison of proportions between groups was completed using {brms} in R and methods described by Andrew Heiss (Heiss, 2023). Data is reported as the mean credible value of the difference between proportions as well as 95% HDI and the probability that the difference in proportions is greater than 0. 
