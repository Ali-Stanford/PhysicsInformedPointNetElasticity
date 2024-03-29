# Physics-Informed PointNet (PIPN) for Linear Elasticity
![pic](./PIPN_elasticity.png) <br>

![pic](./Figure1-1.png) <br>

**Author:** Ali Kashefi (kashefi@stanford.edu)<br>
**Description:** Implementation of Physics-Informed PointNet (PIPN) for weakly-supervised learning of 2D linear elasticity (plane stress) on multiple sets of irregular geometries <br>
**Version:** 1.0 <br>
      
**Citations** <br>
If you use the code, please cite the following journal papers: <br>

<b>#1</b> **[Physics-informed PointNet: On how many irregular geometries can it solve an inverse problem simultaneously? Application to linear elasticity](https://doi.org/10.1615/JMachLearnModelComput.2023050011)**

    @article{kashefi2023PIPNelasticity,
      title={Physics-informed PointNet: On how many irregular geometries can it solve an inverse problem simultaneously? Application to linear elasticity},
      author={Kashefi, Ali and Guibas, Leonidas J and Mukerji, Tapan},
      journal={Journal of Machine Learning for Modeling and Computing},
      volume={4},
      number={4},
      year={2023},
      publisher={Begel House Inc.}}

<b>#2</b> **[Physics-informed PointNet: A deep learning solver for steady-state incompressible flows and thermal fields on multiple sets of irregular geometries](https://doi.org/10.1016/j.jcp.2022.111510)**

    @article{Kashefi2022PIPN, 
      title = {Physics-informed PointNet: A deep learning solver for steady-state incompressible flows and thermal fields on multiple sets of irregular geometries}, 
      journal = {Journal of Computational Physics}, 
      volume = {468},
      pages = {111510}, 
      year = {2022}, 
      issn = {0021-9991}, 
      author = {Ali Kashefi and Tapan Mukerji}}
      
**Abstract** <br>

Regular physics-informed neural networks (PINNs) predict the solution of partial differential equations using sparse labeled data but only over a single domain. On the other hand, fully supervised learning models are first trained usually over a few thousand domains with known solutions (i.e., labeled data), and then predict the solution over a few hundred unseen domains. Physics-informed PointNet (PIPN) is primarily designed to fill this gap between PINNs (as weakly supervised learning models) and fully supervised learning models. In this article, we demonstrate that PIPN predicts the solution of desired partial differential equations over a few hundred domains simultaneously, while it only uses sparse labeled data. This framework benefits fast geometric designs in the industry when only sparse labeled data are available. Particularly, we show that PIPN predicts the solution of a plane stress problem over more than 500 domains with different geometries, simultaneously. Moreover, we pioneer implementing the concept of remarkable batch size (i.e., the number of geometries fed into PIPN at each sub-epoch) into PIPN. Specifically, we try batch sizes of 7, 14, 19, 38, 76, and 133. Additionally, the effect of the PIPN size, symmetric function in the PIPN architecture, and static and dynamic weights for the component of the sparse labeled data in the loss function are investigated.

**Physics-informed PointNet on Wikipedia** <br>
A general description of physics-informed neural networks (PINNs) and their other versions such as PIPN can be found on the following Wikipedia page:<br>
[Physics-informed PointNet (PIPN) for multiple sets of irregular geometries](https://en.wikipedia.org/wiki/Physics-informed_neural_networks#Physics-informed_PointNet_(PIPN)_for_multiple_sets_of_irregular_geometries)

**Physics-informed PointNet Presentation in Machine Learning + X seminar 2022 at Brown University**<br>
In case you are interested, you might watch the recorded machine learning seminar with the topic of PIPN at Brown University using the following link:<br> 
[Video Presentation of PIPN at Brown University](https://www.dropbox.com/s/oafbjl6xaihotqa/GMT20220325-155140_Recording_2560x1440.mp4?dl=0) <br>
[YouTube Video](https://www.youtube.com/watch?v=faeHARnPSVE)

**Questions?** <br>
If you have any questions or need assistance, please do not hesitate to contact Ali Kashefi (kashefi@stanford.edu) via email.

**About the Author** <br>
Please see the author's website: [Ali Kashefi](https://web.stanford.edu/~kashefi/) 
