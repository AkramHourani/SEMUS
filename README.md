# Deep Learning Approaches for Interference Modelling, Detection, and Mitigation for Improving Spaceborne SAR Performance
This research project aims to investigate the interference topic in modern spaceborne SAR systems by first developing a complete open-source spaceborne-SAR simulator/emulator based on near real SAR scenario due to the lake of freely available simulators. The emulator is an RF-level (physical layer) that emulate and analyse the effects of different RFI sources on the final focused spaceborne SAR image, dubbed SEMUS, SAR Emulator for spaceborne application. SEMUS generates the SAR Phase History Data PHD (unfocused raw data) in a low squint case, considering stripmap collection mode for its popularity. Then range Doppler algorithm (RDA) is utilized for generating SAR focused image (single look complex image). A new proposed empirical method is applied for the azimuth matched filter to replace the analytical method in the RDA steps. Furthermore, the emulator is used to inject noise and different types of interference into the PHD to examine the effect of the interference into the emulated SAR raw data and the focused SAR image. 

## Acknowledgments
This work is fully funded by [Smartsat CRC](https://www.smartsatcrc.com/). We are grateful for their funding and commitment to advancing [The Earth observation technology]. Their support has been instrumental in the success of this project.

SmartSat CRC is a consortium of universities and other research organisations, partnered with industry that has been funded by the Australian Government to develop know-how and technologies in advanced telecommunications and IoT connectivity, intelligent satellite systems and Earth observation next generation data services. 

# SEMUS
SAR Emulator for Spaceborne Applications
Please cite:
N. Hendy, F. Kurnia, T. Kraus, M. Bachmann, M. Martorella, R. Evans, M. Zink, H. M. Fayek, A. Al-Hourani, “SEMUS - An Open-Source RF-Level SAR Emulator for Interference Modelling in Spaceborne Applications”, IEEE Transactions on Aerospace and Electronic Systems. (Under review)
