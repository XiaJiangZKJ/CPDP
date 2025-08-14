## Title: The progress in cross-project defect prediction: A revisit and further thinking

## 1. Folder Introduction

- [`CPDP/datasets`](https://github.com/XiaJiangZKJ/CPDP/tree/main/datasets) This folder stores all datasets used in our experiment.

- [`CPDP/scripts`](https://github.com/XiaJiangZKJ/CPDP/tree/main/scripts) This folder stores all codes of ManualDown, other used CPDP methods and our analyzing approaches.

## 2. Execution commands
In order to make it easier to obtain the results, one can run it according to the following command regulation.

(1) run CPDP.R and util.R to get the results of ManualDown. 

(2) run plot.R to get the box plots and line charts between ManualDown and other CPDP methods.

(3) run test.R to analyze the statistical and practical effect of the difference between ManualDown and other CPDP methods.

## 3. Appendix
The following table reports different CPDP approaches' performance using AUC. Additionally, if the method employs a dataset and uses AUC to reflect its effectiveness, the corresponding cell in the table will be filled with the AUC value; if the method employs a dataset but does not use AUC to reflect effectiveness, the corresponding cell will be filled with "✓".

|               |          | DTB-RF-OptADPT | BiLO-CPDP | DEPT-RF | DH-CNN | TBE-HG | DSSDPP | MASTER | SSE   | DP-CLM<sub>E</sub> | DP-GANPT | ManualDown |
|---------------|----------|----------------|-----------|---------|--------|--------|--------|--------|-------|-----------|----------|------------|
| AEEEM         | EQ       | 0.643          | 0.717     |         |        |        | 0.744  | 0.701  |       |           |          | 0.781      |
|               | JDT      | 0.779          | 0.737     |         |        |        | 0.793  | 0.817  |       |           |          | 0.781      |
|               | LC       | 1.000          | 0.709     |         |        |        | 0.790  |        |       |           |          | 0.655      |
|               | ML       | 1.000          | 0.650     |         |        |        | 0.696  |        |       |           |          | 0.692      |
|               | PDE      | 0.739          | 0.685     |         |        |        | 0.722  | 0.729  |       |           |          | 0.717      |
|               | Lucene   | 1.000          |           |         |        |        | 0.814  | 0.781  |       |           |          | 0.652      |
|               | Mylyn    | 1.000          |           |         |        |        | 0.692  | 0.644  |       |           |          | 0.702      |
| JURECZKO      | Ant 1.3  |                |           |         |        |        |        | 0.848  |       |           |          | 0.817      |
|               | Ant 1.5  |                |           | 0.690   |        |        |        |        |       |           |          | 0.788      |
|               | Ant 1.6  |                |           | 0.680   | 0.691  |        |        |        |       | ✓         | ✓        | 0.839       |
|               | Ant 1.7  | 0.969          | 0.800     | 0.740   |        |        | 0.794  |        | 0.789 |           |          | 0.867      |
|               | Camel 1.2|                |           |         |        |        |        |        |       |           |          |            |
|               | Camel 1.4|                |           | 0.610   | 0.574  | 0.567  |        |        |       | ✓         | ✓        | 0.681       |
|               | Camel 1.6| 0.962          | 0.622     | 0.620   |        | 0.614  | 0.653  | 0.613  | 0.619 |           |          | 0.655      |
|               | Ivy 2.0  | 0.736          | 0.807     | 0.740   |        |        | 0.807  | 0.826  | 0.799 | ✓         | ✓        | 0.837      |
|               | Jedit 3.2|                |           |         |        |        |        |        |       |           |          |            |
|               | Jedit 4.0|                |           |         |        | 0.717  |        |        |       |           |          | 0.766      |
|               | Jedit 4.1|                |           | 0.750   | 0.773  | 0.679  |        | 0.800  |       | ✓         | ✓        | 0.805       |
|               | Jedit 4.2|                |           | 0.760   |        | 0.603  |        |        | 0.811 |           |          | 0.830      |
|               | Jedit 4.3| 0.810          | 0.852     | 0.660   |        |        | 0.826  |        |       |           |          | 0.651      |
|               | log4j 1.0|                |           |         |        |        |        |        |       |           |          |            |
|               | log4j 1.1|                |           |         | 0.812  | 0.729  |        |        |       | ✓         | ✓        | 0.799       |
|               | log4j 1.2| 0.798          | 0.842     | 0.640   |        |        | 0.792  | 0.615  | 0.331 |           |          | 0.739      |
|               | Lucene 2.0|                |           |         | 0.685  |        |        |        |       |           |          | 0.734      |
|               | Lucene 2.2|                |           |         |        | 0.597  |        |        |       | ✓         | ✓        | 0.632      |
|               | Lucene 2.4| 0.727          | 0.711     | 0.630   |        | 0.651  | 0.761  |        | 0.673 |           |          | 0.731      |
|               | Poi 1.5  |                |           |         |        |        |        |        |       |           |          |            |
|               | Poi 2.0  |                |           |         |        | 0.586  |        | 0.690  |       |           |          | 0.669      |
|               | Poi 2.5  |                |           | 0.720   |        |        |        |        |       | ✓         | ✓        | 0.697      |
|               | Poi 3.0  | 0.801          | 0.817     | 0.550   | 0.572  | 0.658  | 0.818  |        | 0.74  | ✓         | ✓        | 0.830      |
|               | Synapse 1.1|                |           |         |        |        |        |        |       | ✓         | ✓        |            |
|               | Synapse 1.2| 0.960          | 0.720     | 0.660   | 0.717  | 0.566  | 0.775  | 0.696  | 0.701 | ✓         | ✓        | 0.765      |
|               | Tomcat 6.0|                |           |         |        |        | 0.819  |        |       |           |          | 0.843      |
|               | Velocity 1.4.0|                |           |         |        |        |        | 0.507  |       |           |          | 0.405      |
|               | Velocity 1.5.0|                |           |         |        |        |        |        |       |           |          |            |
|               | Velocity 1.6.0|                |           | 0.590   |        | 0.697  |        |        |       |           |          | 0.694      |
|               | Velocity 1.6.1| 0.739          | 0.702     |         |        |        | 0.676  |        | 0.684 |           |          | 0.714      |
|               | Xalan 2.4.0|                |           |         |        |        |        | 0.778  |       |           |          | 0.798      |
|               | Xalan 2.5.0|                |           |         | 0.725  | 0.608  |        |        |       | ✓         | ✓        | 0.651      |
|               | Xalan 2.6.0|                |           |         |        | 0.641  |        |        |       |           |          | 0.787      |
|               | Xalan 2.7.0| 0.727          | 0.763     |         |        | 0.738  | 0.642  |        | 0.904 |           |          | 0.803      |
|               | Xerces 1.2.0|                |           |         |        |        |        | 0.510  |       |           |          | 0.466      |
|               | Xerces 1.3.0|                |           | 0.470   | 0.658  |        |        |        |       | ✓         | ✓        | 0.753      |
|               | Xerces 1.4.0|                |           | 0.380   |        |        |        |        |       |           |          | 0.755      |
|               | Xerces 1.4.4| 0.701          | 0.716     |         |        |        | 0.659  |        | 0.638 |           |          | 0.681      |
|               | Prop 6.0 |                |           |         |        |        |        | 0.675  |       |           |          | 0.684      |
|               | tomcat   |                |           |         |        |        |        | 0.780  |       |           |          | 0.804      |
| ReLink        | Apache   | 0.837          | 0.748     |         |        |        | 0.736  | 0.725  |       |           |          | 0.776      |
|               | Safe     | 0.962          | 0.799     |         |        |        | 0.814  | 0.808  |       |           |          | 0.834      |
|               | Zxing    | 0.661          | 0.640     |         |        |        | 0.646  | 0.639  |       |           |          | 0.661      |
| GitHub-Python | corefx   |                |           |         |        |        |        | 0.654  |       |           |          | 0.523      |
|               | django   |                |           |         |        |        |        | 0.685  |       |           |          | 0.433      |
|               | nova     |                |           |         |        |        |        | 0.704  |       |           |          | 0.584      |
| JIRA          | Activemq 5.0.0| 0.664          |           |         |        |        | 0.816  | 0.808  |       |           |          | 0.777      |
|               | Derby 10.5.1.1| 0.657          |           |         |        |        | 0.772  | 0.789  |       |           |          | 0.782      |
|               | Groovy 1.6| 0.758          |           |         |        |        | 0.768  | 0.760  |       |           |          | 0.726      |
|               | Hbase 0.94.0| 0.785          |           |         |        |        | 0.768  | 0.782  |       |           |          | 0.763      |
|               | Hive 0.9.0| 0.634          |           |         |        |        | 0.772  | 0.813  |       |           |          | 0.810      |
|               | Jruby 1.1| 0.914          |           |         |        |        | 0.821  | 0.874  |       |           |          | 0.866      |
|               | Wicket 1.3.0| 0.744          |           |         |        |        | 0.818  | 0.824  |       |           |          | 0.827      |
| Eclipse       | Eclipse 2.0| 0.825          |           |         |        |        | 0.760  |        |       |           |          | 0.791      |
|               | Eclipse 2.1| 0.852          |           |         |        |        | 0.689  |        |       |           |          | 0.733      |
|               | Eclipse 3.0| 0.796          |           |         |        |        | 0.717  |        |       |           |          | 0.775      |
| MDP           | JM1      | 0.565          |           |         |        |        | 0.638  |        |       |           |          | 0.663      |
|               | KC3      | 0.836          |           |         |        |        | 0.695  |        |       |           |          | 0.656      |
|               | MC1      | 0.800          |           |         |        |        | 0.749  |        |       |           |          | 0.682      |
|               | MC2      | 0.750          |           |         |        |        | 0.763  |        |       |           |          | 0.672      |
|               | MW1      | 0.835          |           |         |        |        | 0.781  |        |       |           |          | 0.744      |
|               | PC1      | 0.793          |           |         |        |        | 0.726  |        |       |           |          | 0.754      |
|               | PC2      | 0.951          |           |         |        |        | 0.882  |        |       |           |          | 0.562      |
|               | PC3      | 0.648          |           |         |        |        | 0.760  |        |       |           |          | 0.696      |
|               | PC4      | 0.632          |           |         |        |        | 0.676  |        |       |           |          | 0.650      |
|               | PC5      | 0.661          |           |         |        |        | 0.660  |        |       |           |          | 0.713      |
|               | CM1      | 0.762          |           |         |        |        | 0.750  |        |       |           |          | 0.722      |

## 4. Contact us
Mail: 2023141460236@stu.scu.edu.cn