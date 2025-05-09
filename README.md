# Learning a prior for magnetic resonance spectroscopy

This repository contains code and processing files for the final course project in **Advanced Data Science for Biomedical Engineering**. The project adapts the Learned Proximal Network (LPN) framework for denoising Magnetic Resonance Spectroscopy (MRS) signals using synthetic and real-world data.

### Osprey Processing  
`MRS_Osprey_GE_Gaba.m`: MATLAB code to preprocess GE MRS files from Big GABA  
`MRS_Osprey_Siemens_Gaba.m`: MATLAB code to preprocess Siemens MRS files from Big GABA  
`MRS_Osprey_Philips_Gaba.m`: MATLAB code to preprocess Philips MRS files from Big GABA  

### Stage 1  
`stage1_lpn.pth`: Pre-trained Stage 1 model weights  
`Epoch_Logs`: Model output plots every 5 epochs  
`avg_loss_curve.png`: Training loss curve  
`batch_loss_curve.png`: Training loss curve (Batch-wise)  
`gamma_schedule.png`: Gamma decay schedule  
`training_log.csv`: Training log  

### Stage 2  
`stage2_lpn.pth`: Fine-tuned Stage 2 model weights  
`Epoch_Logs`: Model output plots every 5 epochs  
`train_loss_curve.png`: Training loss curve  
`val_loss_curve.png`: Validation loss curve  
`psnr_curve.png`: PSNR curve  
`training_log.csv`: Training log  

## References
[1] Fang, Z., Buchanan, S., & Sulam, J. (2023, October 22). What’s in a prior? Learned proximal networks for inverse problems. arXiv.org. https://arxiv.org/abs/2310.14344  
[2] Gudmundson, A. T., Davies-Jenkins, C. W., Özdemir, İ., Murali-Manohar, S., Zöllner, H. J., Song, Y., Hupfeld, K. E., Schnitzler, A., Oeltzschner, G., Stark, C. E. L., & Edden, R. a. E. (2023). Application of A1H Brain MRS Benchmark Dataset to Deep Learning for Out-of-Voxel artifacts. bioRxiv (Cold Spring Harbor Laboratory). https://doi.org/10.1101/2023.05.08.539813  
[3] Mikkelsen, M., Barker, P. B., Bhattacharyya, P. K., Brix, M. K., Buur, P. F., Cecil, K. M., Chan, K. L., Chen, D. Y., Craven, A. R., Cuypers, K., Dacko, M., Duncan, N. W., Dydak, U., Edmondson, D. A., Ende, G., Ersland, L., Gao, F., Greenhouse, I., Harris, A. D., . . . Edden, R. A. (2017). Big GABA: Edited MR spectroscopy at 24 research sites. NeuroImage, 159, 32–45. https://doi.org/10.1016/j.neuroimage.2017.07.021  
[4] Oeltzschner, G., Zöllner, H. J., Hui, S. C., Mikkelsen, M., Saleh, M. G., Tapper, S., & Edden, R. A. (2020). Osprey: Open-source processing, reconstruction & estimation of magnetic resonance spectroscopy data. Journal of Neuroscience Methods, 343, 108827. https://doi.org/10.1016/j.jneumeth.2020.108827   
[5] Helge, Oeltzschner, G., Davies-Jenkins, C., LaMaster, J., Craven, A. R., Richardedden, P, V. S., Mikkelsen, M., Helmick, C., Hui, S., Meganaforbes, Raumhein, & Evans, J. (2024). schorschinho/osprey: Osprey v.2.9.0. Zenodo. https://doi.org/10.5281/zenodo.14226197




