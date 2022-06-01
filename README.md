#Installation instructions / Run Instructions

I used the MASS, survival, randomForestSRC, alabama, and survminer package for R and you'll notice code source "ParametricEst_function.R", "cox survival function.R", "stacking model.R", "Covid model.R", "ISSE.R", and "ROC.R".

Step 1: Install all package in your R (MASS, survival, randomForestSRC,alabama, survminer).

Step 2: Run "ParametricEst_function.R", "cox survival function.R", "code stacking model.R" .

Step 2a: Open "Covid model.R" in R and start inputting your data (covid dataset) by dat1=read.delim("clipboard")

Step 2b: Run all code in "Covid model.R" after inputing data code.

Step 3: Run all code "ROC.R" and then change the marker by type model.

Step 4: write summary(mod.ln) and summary(mod.cox) for result of each model.

Step 5: write stacked.est$alphas to weight from each model for stacked model.

#Requirements

This was created and tested on an R X64. 4.0.3.
