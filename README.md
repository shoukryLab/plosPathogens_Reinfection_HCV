# plosPathogens_Reinfection_HCV
For information on the projects pertinent to the code uploaded here, please refer to the following publication: Differential Immune Transcriptomic Profiles between Vaccinated and Resolved HCV Reinfected Subjects'. In short, the aformentioned paper looks into the transcriptional profile of subjects re-infected with HCV in a cohort of PWID. 

A previous publication in collaboration with the same group looks at a similar mechanism, however with subjects primary infection to HCV, titled: 'Longitudinal transcriptomic characterization of the immune response to acute hepatitis C virus infection in patients with spontaneous viral clearance'. Expression profiles from both of these publications were taken into consideration in order to dig deeper into the underlying mechanisms in play during primary and secondary HCV infection; for more information, you can refer to both of the mentioned publications.

The eset containing the metadata, as well as count matrix for the aforementioned studies has been uploaded in order to facilitate the analysis steps for those who wish to reproduce our results. 

In the metadata of the eset - which can be obtained by running pData(esetRaw), after loading the data - most of the column names are self explanatory. The one used in the linear model - Class, less so - however is paramount to our analysis, and consists of a concatenation of two factors: `Time.point` and `Viral.load`, separated by an underscore -> e.g. "LAI_neg", represents a Late Acute I timepoint (i.e. Late Acute primary infection), who is NO LONGER INFECTED with HCV.



