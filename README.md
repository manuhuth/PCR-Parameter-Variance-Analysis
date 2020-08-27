# Student Project: Analysis of the Variance of the Principal Component Regression Coefficients and the Estimated Outcome
### Does Knowing the True Variance Covariance Matrix Decrease the Variance?

## This Repository
This repository contains my approach to find evidence how the variance of the principal component regression coefficients and the estimated outcome changes, if the true variance covariance matrix is unknown. I have created an extensive data generating process to approximate real world data and derived theory to set up a proper simulation study.
The main parts of the analysis are the created R functions inside the folder [R](https://github.com/manuhuth/PCR-Parameter-Variance-Analysis/tree/master/R), the [notebook](https://github.com/manuhuth/PCR-Parameter-Variance-Analysis/blob/master/Notebook.ipynb) that is mainly created to reproduce the conducted simulation study and thus does not include proofs of theorems used to derive the set up of the simulation study or the data generating process and the [pdf](https://github.com/manuhuth/PCR-Parameter-Variance-Analysis/blob/master/PDF.pdf) hat is similar to the notebook but contains additionally proofs of theorems used to derive the set up of the simulation study or the data generating process. In the notebook it is referred to the page in the paper, when technical details are left out to increase readability. 

The ensure that every image or format is displayed properly, I recommend to download this notebook from its repository on [GitHub](https://github.com/manuhuth/PCR-Parameter-Variance-Analysis). However, the following badges allow easy access to the project's notebook

[![nbviewer](https://camo.githubusercontent.com/bfeb5472ee3df9b7c63ea3b260dc0c679be90b97/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f72656e6465722d6e627669657765722d6f72616e67652e7376673f636f6c6f72423d66333736323626636f6c6f72413d346434643464)](https://nbviewer.jupyter.org/github/manuhuth/PCR-Parameter-Variance-Analysis/blob/master/Notebook.ipynb)

<a href="https://mybinder.org/v2/gh/manuhuth/PCR-Parameter-Variance-Analysis/master?filepath=notebook"
    target="_parent">
    <img align="center"
       src="https://mybinder.org/badge_logo.svg"
       width="109" height="20">
</a>



## Continous integration with Travis CI to ensure Reproducibility
To ensure reproducibility, I have integrated Travis CI. The build history can be found here [![Build Status](https://travis-ci.org/HumanCapitalAnalysis/microeconometrics-course-project-manuhuth.svg?branch=master)](https://travis-ci.org/github/manuhuth/PCR-Parameter-Variance-Analysis)


## References

> Bailey, M. J. (2010). " Momma's got the pill": how Anthony Comstock and Griswold v. Connecticut shaped US childbearing. American economic review, 100(1), 98-129. [Link](https://www.aeaweb.org/articles?id=10.1257/aer.100.1.98)

> Björklund, A., & Kjellström, C. (2002). Estimating the return to investments in education: how useful is the standard Mincer equation?. Economics of Education Review, 21(3), 195-210. [Link](https://www.sciencedirect.com/science/article/abs/pii/S0272775701000036)

> Blundell, R., Dearden, L., & Sianesi, B. (2005). Evaluating the effect of education on earnings: models, methods and results from the National Child Development Survey. Journal of the Royal Statistical Society: Series A (Statistics in Society), 168(3), 473-512. [Link](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-985X.2004.00360.x?casa_token=8XTSrhMvcoIAAAAA%3Ar0ZECHQIWsbtdynj4kZZ_R_-HSDkUKPlkLvS8GF9whkNF584aPmn6nHGR4cZXOOZTVLQQu_-9E8VunWZ)

> Blundell, R., Costa Dias, M., Meghir, C., & Shaw, J. (2016). Female labor supply, human capital, and welfare reform. Econometrica, 84(5), 1705-1753. [Link](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA11576?casa_token=T7_8JfSc3V0AAAAA:Rq3dQoE4CwRQ5bmOvmvwz1RfUnzp7wYh3LRtrMZsrJjKhrefbBItL0gdFcdZLiYVB-33AUya90S8wTfy)

> Blundell, R., Pistaferri, L., & Saporta-Eksten, I. (2016). Consumption inequality and family labor supply. American Economic Review, 106(2), 387-435. [Link](https://www.aeaweb.org/articles?id=10.1257/aer.20121549)

> British Government. School Starting Age. visited at 10.08.2020 at 19.07. [Link](https://www.gov.uk/schools-admissions/school-starting-age)

> British Government. School Leaving Age. visited at 10.08.2020 at 19.17. [Link](https://www.gov.uk/know-when-you-can-leave-school)

> British Government. State Pension Age. visited at 10.08.2020 at 18.13. [Link](https://www.gov.uk/state-pension-age)

> Card, D. (2001). Estimating the return to schooling: Progress on some persistent econometric problems. Econometrica, 69(5), 1127-1160. [Link](https://onlinelibrary.wiley.com/doi/full/10.1111/1468-0262.00237?casa_token=Fb2oj5pOnrsAAAAA%3AQmDyCEKfwLAsw9z4b1JZgjqWnCgpYw0n49-ljsiGMyx1va5NnTVvAdUL_U907RhsL_EtJ10m3V2iDdfV)

> Carneiro, P., Heckman, J. J., & Vytlacil, E. J. (2011). Estimating marginal returns to education. American Economic Review, 101(6), 2754-81. [Link](https://www.aeaweb.org/articles?id=10.1257/aer.101.6.2754)

> Cygan-Rehm, K., & Maeder, M. (2013). The effect of education on fertility: Evidence from a compulsory schooling reform. Labour Economics, 25, 35-48. [Link](https://www.sciencedirect.com/science/article/abs/pii/S0927537113000584)

> Consul, P. C., & Jain, G. C. (1973). A generalization of the Poisson distribution. Technometrics, 15(4), 791-799. [Link](https://www.tandfonline.com/doi/abs/10.1080/00401706.1973.10489112)

> Davis-Kean, P. E. (2005). The influence of parent education and family income on child achievement: the indirect role of parental expectations and the home environment. Journal of family psychology, 19(2), 294. [Link](https://psycnet.apa.org/buy/2005-06518-016)

> Fan, X., Seshadri, A., & Taber, C. (2015). Estimation of a life-cycle model with human capital, labor supply and retirement. Australia. University of New South Wales. [Link](https://www.ssc.wisc.edu/~aseshadr/WorkingPapers/FST.pdf)

> Friedman, J., Hastie, T., & Tibshirani, R. (2001). The elements of statistical learning (Vol. 1, No. 10). New York: Springer series in statistics. [Link](https://psycnet.apa.org/buy/2005-06518-016)

> Hansen, K. T., Heckman, J. J., & Mullen, K. J. (2004). The effect of schooling and ability on achievement test scores. Journal of econometrics, 121(1-2), 39-98. [Link](https://www.sciencedirect.com/science/article/abs/pii/S0304407603002598)

> Heckman, J. J., Lochner, L. J., & Todd, P. E. (2006). Earnings functions, rates of return and treatment effects: The Mincer equation and beyond. Handbook of the Economics of Education, 1, 307-458. [Link](https://www.sciencedirect.com/science/article/pii/S1574069206010075)

> Holmlund, B., Liu, Q., & Nordström Skans, O. (2008). Mind the gap? Estimating the effects of postponing higher education. Oxford Economic Papers, 60(4), 683-710. [Link](https://academic.oup.com/oep/article-abstract/60/4/683/2362081)

> Human Development Report (2020). United Nations Development Programme. [Link](http://hdr.undp.org/en/indicators/103006)

> James, G., Witten, D., Hastie, T., & Tibshirani, R. (2013). An introduction to statistical learning (Vol. 112, p. 18). New York: springer. [Link](https://link.springer.com/book/10.1007%2F978-1-4614-7138-7)

> Jolliffe, I. T. (1986). Principal components in regression analysis. In Principal component analysis (pp. 129-155). Springer, New York, NY. [Link](https://link.springer.com/chapter/10.1007/978-1-4757-1904-8_8)

> Jung, R. C., & Tremayne, A. R. (2011). Convolution‐closed models for count time series with applications. Journal of Time Series Analysis, 32(3), 268-280. [Link](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9892.2010.00697.x?casa_token=Z_aij9JL6U0AAAAA%3APnfZIk7h3P8Jzs4pVvW-Hmd1mswaOAMW-5KbW6G6YV3uEQiXvbrT7ZH8pGHJqiSed7ofNfkwSbgb2hfx)

> Lemieux, T. (2006). The “Mincer equation” thirty years after schooling, experience, and earnings. In Jacob Mincer a pioneer of modern labor economics (pp. 127-145). Springer, Boston, MA. [Link](https://link.springer.com/chapter/10.1007/0-387-29175-X_11)

> Mincer, J. (1974). Schooling, Experience, and Earnings. Human Behavior & Social Institutions No. 2. [Link](https://eric.ed.gov/?id=ED103621)

> Ku, W., Storer, R. H., & Georgakis, C. (1995). Disturbance detection and isolation by dynamic principal component analysis. Chemometrics and intelligent laboratory systems, 30(1), 179-196. [Link](https://www.sciencedirect.com/science/article/pii/0169743995000763)

> Office for National Statistics. Population Estimates. visited at 10.08.2020 at 19.02. [Link](https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland/mid2015/ukandregionalpopulationestimates18382015.zip)

> Shlens, J. (2014). A tutorial on principal component analysis. arXiv preprint arXiv:1404.1100. [Link](https://www.cs.cmu.edu/~elaw/papers/pca.pdf)

> Willis, R. J. (1973). A new approach to the economic theory of fertility behavior. Journal of political Economy, 81(2, Part 2), S14-S64. [Link](https://www.journals.uchicago.edu/doi/abs/10.1086/260152?journalCode=jpe)
