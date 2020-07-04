# Introduction

The **extended generalized classic yield function** (XGC) is convex and C2 continuous criterion. It can describe basic shear strength, tensile strength cut-off, compressive strength cap, and the impact of intermediate principal stress.

<center> 
    <img src="Image/FigOverall.jpg" width="300"> <br> 
    <div style="display: inline-block;">
    Fig.1 the XGC yield surface in principal stress space</div> 
</center>


# Parameters

There are five parameters in XGC.

<img src="svgs/9c8f92472352eba54f4bbf577321e74a.svg?invert_in_darkmode" align=middle width=52.7076pt height=24.65759999999998pt/> are basic shear strength parameters. Their physical meanings are the cohesion and friction angle at the Lode angle <img src="svgs/fd1203f97cc7ecee568cec073736dcef.svg?invert_in_darkmode" align=middle width=38.31036pt height=22.831379999999992pt/>.

<img src="svgs/ae8b681b4ade976a90acb5cd6297df2a.svg?invert_in_darkmode" align=middle width=63.133455pt height=24.65759999999998pt/> is a parameter for the impact of the intermediate principal stress. A larger <img src="svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode" align=middle width=10.165650000000005pt height=22.831379999999992pt/> indicates a smaller impact, and <img src="svgs/db4c2f097a22313356725fc8d816eb43.svg?invert_in_darkmode" align=middle width=40.302405pt height=22.831379999999992pt/> recovers to the MC criterion, where the shear strength only controlled by only the major and the minor principal stresses.
    
<img src="svgs/da3c8f0f216a8a53e5dacd2c52d301dd.svg?invert_in_darkmode" align=middle width=117.18200999999998pt height=21.18732pt/> represents the tensile strength. Its physical meaning is the maximum mean stress at the condition of zero deviatoric stress. (Note: tensile strength is always suggested, because it can avoid the singularities at <img src="svgs/568946085eeaf880fcbe810adcf7af1a.svg?invert_in_darkmode" align=middle width=40.83321pt height=22.46574pt/>)
    
<img src="svgs/9b9415527ca0c0ac5204af6bf2f0eb09.svg?invert_in_darkmode" align=middle width=48.222899999999996pt height=21.18732pt/> (or not valid) represents the compressive strength cap. Its physical meaning is the maximum negative mean stress at the condition of zero deviatoric stress. 
    
The yield function can also provide a close approximation to **MC criterion**, to avoid the issue of discontinuities gradients. In this case, the suggested parameters are <img src="svgs/1b207514b624bf80a53ab868a5fa0bde.svg?invert_in_darkmode" align=middle width=351.11785499999996pt height=22.831379999999992pt/>. However, for most cases, a tensile strength <img src="svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.936155500000004pt height=20.222069999999988pt/> close to zero is more suitable for soils.

# Calibration
<img src="svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.936155500000004pt height=20.222069999999988pt/> and <img src="svgs/0e51a2dede42189d77627c4d742822c3.svg?invert_in_darkmode" align=middle width=14.433210000000003pt height=14.155350000000013pt/> are determined by their physical meanings.

<img src="svgs/154cee6b660d527c4cda27f69e04b4a6.svg?invert_in_darkmode" align=middle width=70.178955pt height=24.65759999999998pt/> can be calibrated by 
triaxial compression and triaxial extension tests. 

- <img src="svgs/70ea4c8e89ff931f610ee360d0db11b8.svg?invert_in_darkmode" align=middle width=30.310335000000002pt height=14.155350000000013pt/>: the cohesion measured by triaxial compression test;
- <img src="svgs/7f079413ee7859f4236782ef401e03bc.svg?invert_in_darkmode" align=middle width=33.949905pt height=14.155350000000013pt/>: the friction angle measured by triaxial compression test;
- <img src="svgs/a78f161e5f0443b60d6708691f4aca2b.svg?invert_in_darkmode" align=middle width=23.504745000000003pt height=14.155350000000013pt/>: the friction angle measured by triaxial extension test;

The algorithm is as below.

<center>
    <img src="Image/TableCalibration.jpg" width ="400">
</center>

# Sample codes

[sub_EGC.m]: sub_EGC.m

The programming of yield value, gradient, and Hessian is provided in '[sub_EGC.m]' 