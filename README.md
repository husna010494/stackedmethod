#Installation instructions / Run Instructions

I used the MASS, survival, randomForestSRC, alabama, and survminer package for R and you'll notice code source "Parametric survival.R", "cox survival function.R", "code stacking model.R", "ISSE.R", and "Kaplan Meier curve.R" in .rar files.

Step 1: Install all package in your R (MASS, survival, randomForestSRC,alabama, survminer).

Step 2: Run "Parametric survival.R", "cox survival function.R", "code stacking model.R" .

Step 2a: Open "code stacking model.R" in R and start inputting your data (covid2 (1)) for .csv files in R (I already wrote the code in code stacking model.R files)

Step 2b: Run all code in "code stacking model.R" after inputing data code.

Step 3: Run all code "ISSE.R".

Step 4: write summary(mod.ln) and summary(mod.cox) for result of each model.

step 5: write "error1" for ISSE log-normal survival model, write "error2" for ISSE cox survival model, write "error3" for ISSE random forest survival model,
and write "error4" for ISSE stacked model.

step 6: run all code "Kaplan Meier curve.R" to get the Kaplan-Meier curve.

#Requirements

This was created and tested on an R X64. 4.0.3.
