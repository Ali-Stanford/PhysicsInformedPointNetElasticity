# Physics-Informed PointNet (PIPN) for Linear Elasticity

**Author:** Ali Kashefi (kashefi@stanford.edu)<br>
**Description:** Implementation of physics-informed PointNet (PIPN) for weakly-supervised learning of 2D linear elasticity (plane stress) on multiple sets of irregular geometries <br>
**Version:** 1.0 <br>

**Citation** <br>
If you use the code, plesae cite the following journal papers: <br>

1. **[Physics-informed PointNet: On how many irregular geometries can it solve an inverse problem simultaneously? Application to linear elasticity](https://arxiv.org/pdf/2303.13634.pdf)**

@article{kashefi2023PIPNelasticity, <br>
title={Physics-informed PointNet: On how many irregular geometries can it solve an inverse problem simultaneously? Application to linear elasticity}, <br>
  author={Kashefi, Ali and Mukerji, Tapan}, <br>
  journal={arXiv preprint arXiv:2303.13634}, <br>
  year={2023}} <br>

2. **[Physics-informed PointNet: A deep learning solver for steady-state incompressible flows and thermal fields on multiple sets of irregular geometries](https://doi.org/10.1016/j.jcp.2022.111510)**

@article{Kashefi2022PIPN, <br>
title = {Physics-informed PointNet: A deep learning solver for steady-state incompressible flows and thermal fields on multiple sets of irregular geometries}, <br>
journal = {Journal of Computational Physics}, <br>
volume = {468}, <br>
pages = {111510}, <br>
year = {2022}, <br>
issn = {0021-9991}, <br>
author = {Ali Kashefi and Tapan Mukerji}} <br>

**Abstract** <br>

Regular physics-informed neural networks (PINNs) predict the solution of partial differential equations using sparse labeled data but only over a single domain. On the other hand, fully supervised learning models are first trained usually over a few thousand domains with known solutions (i.e., labeled data) and then predict the solution over a few hundred unseen domains. Physics-informed PointNet (PIPN) is primarily designed to fill this gap between PINNs (as weakly supervised learning models) and fully supervised learning models. In this article, we demonstrate that PIPN predicts the solution of desired partial differential equations over a few hundred domains simultaneously, while it only uses sparse labeled data. This framework benefits fast geometric designs in the industry when only sparse labeled data are available. Particularly, we show that PIPN predicts the solution of a plane stress problem over more than 500 domains with different geometries, simultaneously. Moreover, we pioneer implementing the concept of remarkable batch size (i.e., the number of geometries fed into PIPN at each sub-epoch) into PIPN. Specifically, we try batch sizes of 7, 14, 19, 38, 76, and 133. Additionally, the effect of the PIPN size, symmetric function in the PIPN architecture, and static and dynamic weights for the component of the sparse labeled data in the loss function are investigated.

**Physics-informed PointNet on Wikipedia** <br>
A general description of physics-informed neural networks (PINNs) and its other versions such as PIPN can be found in the following Wikipedia page:<br>
[Physics-informed PointNet (PIPN) for multiple sets of irregular geometries](https://en.wikipedia.org/wiki/Physics-informed_neural_networks#Physics-informed_PointNet_(PIPN)_for_multiple_sets_of_irregular_geometries)

**Physics-informed PointNet Presentation in Machine Learning + X seminars 2022 at Brown University**<br>
In case of your interest, you might watch the recorded machine learning seminar with the topic of PIPN at Brown University using the following link:<br> 
[Video Presentation of PIPN at Brown University](https://www.dropbox.com/s/oafbjl6xaihotqa/GMT20220325-155140_Recording_2560x1440.mp4?dl=0)


**Questions?** <br>
If you have any questions or need assistance, please do not hesitate to contact Ali Kashefi (kashefi@stanford.edu) via email. 
