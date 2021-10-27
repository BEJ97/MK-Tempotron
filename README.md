# MK-Tempotron
Efficient Motion Symbol Detection and Multi-kernel Learning for Event-based Object Recognition(TCDS 2021)

![imgage](https://github.com/BEJ97/MK-Tempotron/blob/main/image/framework%20of%20MK-Tempotron.png)
## Summary
This is the implemtation code for the following paper.

cite as:

Y. Chen, W. Bai, Q. Huang and J. Xiao, "Efficient Motion Symbol Detection and Multi-kernel Learning for AER Object Recognition," in IEEE Transactions on Cognitive and Developmental Systems, doi: 10.1109/TCDS.2021.3122131.

Bibtex:

```
@ARTICLE{9585046,
author={Chen, Yunhua and Bai, Weijie and Huang, Qingkun and Xiao, Jinsheng},
journal={IEEE Transactions on Cognitive and Developmental Systems},
title={Efficient Motion Symbol Detection and Multi-kernel Learning for AER Object Recognition},
year={2021},
volume={},
number={},
pages={1-1},
doi={10.1109/TCDS.2021.3122131}}
```

## MK-Tempotron

[Link to paper](https://ieeexplore.ieee.org/document/9585046)

Address event representation (AER) vision sensors process visual information by monitoring changes in light intensity and generating event streams. In this paper, a threshold mechanism based on the statistical results of the peak membrane potentials, is proposed to improve the accuracy and efficiency of the motion symbol detection (MSD) method based on the leaky integrate-and-fire (LIF) neuron model and a peak spiking monitoring unit (PSMU). And a multi-kernel learning algorithm based on the tempotron rule, namely MK-tempotron, is proposed to improve the antinoise performance of the classifier. In MK-tempotron, different kernels are applied to calculate the post-synaptic membrane potentials (PSPs) according to different input spiking patterns, to counteract the impact of noise on the spiking activities, so that the synaptic weights can be changed in the direction of correct firing. We verified the effectiveness of the threshold mechanism and multi-kernel learning on a variety of data sets. The experimental results show that among the several recent algorithms based on temporal encoding and SNN classifier, which have advantages in power consumption and network latency, our method achieved the best accuracy on MNIST-DVS (100ms) and N-CARS, and it also achieved competitive results on N-MNIST, CIFAR10-DVS, Posture-DVS, Poker-DVS and MNIST-DVS (200ms/500ms/2000ms).

## Code Implementation
### Requirement
```
matlab
```
### Preparations:
1. Download Datasets 
2. Transfer raw datasets to .mat files (or download our processed mat files following the instruction in readme.txt in code folder)
### Running examples:
Change the Matlab directory to the code folder and execute the Matlab script Top_on***.m based on which datasets you want to test.
