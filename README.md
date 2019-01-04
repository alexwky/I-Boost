# I-Boost
I-Boost is a statistical boosting method that integrates multiple types of high-dimensional genomics data with clinical data for predicting survival time. It is implemented in the R-package **IBoost**.

**IBoost** relies on the R-package **glmnet**, which is hosted on CRAN.
```
install.packages("glmnet")
```

**IBoost** can be installed from github directly:
```
install.packages("devtools")
library(devtools)
install_github("alexwky/I-Boost")
```
**IBoost()** is the main function that implements the I-Boost procedure. Details about the function can be found in the user manual:
```
?IBoost
```

## Contact
**Kin Yau Wong** <kin-yau.wong@polyu.edu.hk>
