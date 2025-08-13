## Title: The progress in cross-project defect prediction: A revisit and further thinking

## 1. Folders Introduction

- [`CPDP/datasets`]() This folder stores all datasets used in our experiment.

- [`CPDP/scripts`]() This folder stores all codes of ManualDown, other used CPDP methods and our analyzing approaches.

## 2. Execution commands
In order to make it easier to obtain the results, one can run it according to the following command regulation.
### (1) run CPDP.R and util.R to get the results of ManualDown. 
### (2) run plot.R to get the box plots and line charts between ManualDown and other CPDP methods.
### (3) run test.R to analyze the statistical and practical effect of the difference between ManualDown and other CPDP methods.

## 3. Appendix
The following table reports different CPDP approaches' performance using AUC. Additionally, if the method employs a dataset and uses AUC to reflect its effectiveness, the corresponding cell in the table will be filled with the AUC value; if the method employs a dataset but does not use AUC to reflect effectiveness, the corresponding cell will be filled with "✓".

|               |          | DTB-RF-OptADPT | BiLO-CPDP | DEPT-RF | DH-CNN | TBE-HG | DSSDPP | MASTER | SSE   | DP-CLM<sub>E</sub> | DP-GANPT |
|---------------|----------|----------------|-----------|---------|--------|--------|--------|--------|-------|-----------|----------|
| AEEEM         | EQ       | 0.643          | 0.717     |         |        |        | 0.744  | 0.701  |       |           |          |
|               | JDT      | 0.779          | 0.737     |         |        |        | 0.793  | 0.817  |       |           |          |
|               | LC       | 1.000          | 0.709     |         |        |        | 0.790  |        |       |           |          |
|               | ML       | 1.000          | 0.650     |         |        |        | 0.696  |        |       |           |          |
|               | PDE      | 0.739          | 0.685     |         |        |        | 0.722  | 0.729  |       |           |          |
|               | Lucene   | 1.000          |           |         |        |        | 0.814  | 0.781  |       |           |          |
|               | Mylyn    | 1.000          |           |         |        |        | 0.692  | 0.644  |       |           |          |
| JURECZKO      | Ant 1.3  |                |           |         |        |        |        | 0.848  |       |           |          |
|               | Ant 1.5  |                |           | 0.690   |        |        |        |        |       |           |          |
|               | Ant 1.6  |                |           | 0.680   | 0.691  |        |        |        |       | ✓         | ✓        |
|               | Ant 1.7  | 0.969          | 0.800     | 0.740   |        |        | 0.794  |        | 0.789 |           |          |
|               | Camel 1.2|                |           |         |        |        |        |        |       |           |          |
|               | Camel 1.4|                |           | 0.610   | 0.574  | 0.567  |        |        |       | ✓         | ✓        |
|               | Camel 1.6| 0.962          | 0.622     | 0.620   |        | 0.614  | 0.653  | 0.613  | 0.619 |           |          |
|               | Ivy 2.0  | 0.736          | 0.807     | 0.740   |        |        | 0.807  | 0.826  | 0.799 | ✓         | ✓        |
|               | Jedit 3.2|                |           |         |        |        |        |        |       |           |          |
|               | Jedit 4.0|                |           |         |        | 0.717  |        |        |       |           |          |
|               | Jedit 4.1|                |           | 0.750   | 0.773  | 0.679  |        | 0.800  |       | ✓         | ✓        |
|               | Jedit 4.2|                |           | 0.760   |        | 0.603  |        |        | 0.811 |           |          |
|               | Jedit 4.3| 0.810          | 0.852     | 0.660   |        |        | 0.826  |        |       |           |          |
|               | log4j 1.0|                |           |         |        |        |        |        |       |           |          |
|               | log4j 1.1|                |           |         | 0.812  | 0.729  |        |        |       | ✓         | ✓        |
|               | log4j 1.2| 0.798          | 0.842     | 0.640   |        |        | 0.792  | 0.615  | 0.331 |           |          |
|               | Lucene 2.0|                |           |         | 0.685  |        |        |        |       |           |          |
|               | Lucene 2.2|                |           |         |        | 0.597  |        |        |       | ✓         | ✓        |
|               | Lucene 2.4| 0.727          | 0.711     | 0.630   |        | 0.651  | 0.761  |        | 0.673 |           |          |
|               | Poi 1.5  |                |           |         |        |        |        |        |       |           |          |
|               | Poi 2.0  |                |           |         |        | 0.586  |        | 0.690  |       |           |          |
|               | Poi 2.5  |                |           | 0.720   |        |        |        |        |       | ✓         | ✓        |
|               | Poi 3.0  | 0.801          | 0.817     | 0.550   | 0.572  | 0.658  | 0.818  |        | 0.74  | ✓         | ✓        |
|               | Synapse 1.1|                |           |         |        |        |        |        |       | ✓         | ✓        |
|               | Synapse 1.2| 0.960          | 0.720     | 0.660   | 0.717  | 0.566  | 0.775  | 0.696  | 0.701 | ✓         | ✓        |
|               | Tomcat 6.0|                |           |         |        |        | 0.819  |        |       |           |          |
|               | Velocity 1.4.0|                |           |         |        |        |        | 0.507  |       |           |          |
|               | Velocity 1.5.0|                |           |         |        |        |        |        |       |           |          |
|               | Velocity 1.6.0|                |           | 0.590   |        | 0.697  |        |        |       |           |          |
|               | Velocity 1.6.1| 0.739          | 0.702     |         |        |        | 0.676  |        | 0.684 |           |          |
|               | Xalan 2.4.0|                |           |         |        |        |        | 0.778  |       |           |          |
|               | Xalan 2.5.0|                |           |         | 0.725  | 0.608  |        |        |       | ✓         | ✓        |
|               | Xalan 2.6.0|                |           |         |        | 0.641  |        |        |       |           |          |
|               | Xalan 2.7.0| 0.727          | 0.763     |         |        | 0.738  | 0.642  |        | 0.904 |           |          |
|               | Xerces 1.2.0|                |           |         |        |        |        | 0.510  |       |           |          |
|               | Xerces 1.3.0|                |           | 0.470   | 0.658  |        |        |        |       | ✓         | ✓        |
|               | Xerces 1.4.0|                |           | 0.380   |        |        |        |        |       |           |          |
|               | Xerces 1.4.4| 0.701          | 0.716     |         |        |        | 0.659  |        | 0.638 |           |          |
|               | Prop 6.0 |                |           |         |        |        |        | 0.675  |       |           |          |
|               | tomcat   |                |           |         |        |        |        | 0.780  |       |           |          |
| ReLink        | Apache   | 0.837          | 0.748     |         |        |        | 0.736  | 0.725  |       |           |          |
|               | Safe     | 0.962          | 0.799     |         |        |        | 0.814  | 0.808  |       |           |          |
|               | Zxing    | 0.661          | 0.640     |         |        |        | 0.646  | 0.639  |       |           |          |
| GitHub-Python | corefx   |                |           |         |        |        |        | 0.654  |       |           |          |
|               | django   |                |           |         |        |        |        | 0.685  |       |           |          |
|               | nova     |                |           |         |        |        |        | 0.704  |       |           |          |
| JIRA          | Activemq 5.0.0| 0.664          |           |         |        |        | 0.816  | 0.808  |       |           |          |
|               | Derby 10.5.1.1| 0.657          |           |         |        |        | 0.772  | 0.789  |       |           |          |
|               | Groovy 1.6| 0.758          |           |         |        |        | 0.768  | 0.760  |       |           |          |
|               | Hbase 0.94.0| 0.785          |           |         |        |        | 0.768  | 0.782  |       |           |          |
|               | Hive 0.9.0| 0.634          |           |         |        |        | 0.772  | 0.813  |       |           |          |
|               | Jruby 1.1| 0.914          |           |         |        |        | 0.821  | 0.874  |       |           |          |
|               | Wicket 1.3.0| 0.744          |           |         |        |        | 0.818  | 0.824  |       |           |          |
| Eclipse       | Eclipse 2.0| 0.825          |           |         |        |        | 0.760  |        |       |           |          |
|               | Eclipse 2.1| 0.852          |           |         |        |        | 0.689  |        |       |           |          |
|               | Eclipse 3.0| 0.796          |           |         |        |        | 0.717  |        |       |           |          |
| MDP           | JM1      | 0.565          |           |         |        |        | 0.638  |        |       |           |          |
|               | KC3      | 0.836          |           |         |        |        | 0.695  |        |       |           |          |
|               | MC1      | 0.800          |           |         |        |        | 0.749  |        |       |           |          |
|               | MC2      | 0.750          |           |         |        |        | 0.763  |        |       |           |          |
|               | MW1      | 0.835          |           |         |        |        | 0.781  |        |       |           |          |
|               | PC1      | 0.793          |           |         |        |        | 0.726  |        |       |           |          |
|               | PC2      | 0.951          |           |         |        |        | 0.882  |        |       |           |          |
|               | PC3      | 0.648          |           |         |        |        | 0.760  |        |       |           |          |
|               | PC4      | 0.632          |           |         |        |        | 0.676  |        |       |           |          |
|               | PC5      | 0.661          |           |         |        |        | 0.660  |        |       |           |          |
|               | CM1      | 0.762          |           |         |        |        | 0.750  |        |       |           |          |