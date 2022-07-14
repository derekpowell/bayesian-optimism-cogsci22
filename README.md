# A descriptive Bayesian account of optimism in belief revision

## Materials

Life event items were as follows:

| item                                                     | true probability (%) |
|----------------------------------------------------------|-----------------|
| be a victim of a car break-in                            | 76              |
| die before the age of 80                                 | 41              |
| suffer heart failure                                     | 31              |
| die before the age of 70                                 | 18              |
| suffer a stroke                                          | 17              |
| be physically assaulted by a stranger                    | 26              |
| be diagnosed with any kind of cancer                     | 40              |
| become divorced                                          | 45              |
| have cataracts                                           | 65              |
| be diagnosed with skin cancer                            | 3               |
| become a homeowner                                       | 82              |
| live past 90 years old                                   | 11              |
| maintain a healthy weight                                | 58              |
| have an achievement recognized by the press              | 3               |
| have the financial security to retire comfortably at 65  | 22              |
| live past 100 years old                                  | 1               |
| be in a marriage that lasts at least 10 years            | 61              |
| be earning more than $75k/yr within 5 years              | 16              |
| be earning more than $100k/yr within 5 years             | 9               |
| make it through all of next year without catching a cold | 23              |

## Reproducibility notes

To reproduce the manuscript and all analyses, follow the following steps after cloning this repository.

**RECOMMENDED**: Install [Docker](https://www.docker.com/) and utilize `cogdatasci/rstudio` docker container. Recommend you run with the following options:
```bash
docker run -d -p 8787:8787 -v "`pwd`":/home/rstudio/working \
 -e PASSWORD=my_password_here cogdatasci/rstudio
 ```
then navigate to `localhost:8787` in your browser to access the rstudio interface.
