# EME 和 3D-FDTD 仿真均匀布拉格光栅时的对比



## 背景

我们常需要准确的仿真均匀的短长度的光栅。

因为光栅许多特性，仅由其单元周期光栅的结构所决定，如$\kappa$, 不同耦合模式对应的反射带中心，等等。

- 仿真短长度的光栅，并通过其反射谱的反射率峰值，计算出光栅的 $\kappa$， 通过以下公式 (Eq. 4.32 in Silicon Photonics Design,  Lukas Chrostowski)

$$
R_{peak} = tanh^2(\kappa L)
$$

- 仿真短长度的光栅，可先预估某个反射谱中心。
- 仿真短长度的光栅，可先预估其不同耦合模式的反射带的位置 （如模式转换的光栅辅助反向耦合器）。

当我知道了具有某个单元光栅周期结构的$\kappa$, 我们便可使用基于CMT的传输矩阵法，进一步设计我们光栅的光谱。



## EME 和 3D-FDTD仿真结果对比的总结

光栅反射谱的计算方法，EME 和 3D-FDTD。两种方法均可得到较为准确的结果。

但是EME，可通过计算**单个周期单元的光栅**的模式 ，便可实现整体器件的光谱的计算，**更加快速** （**比3D-FDTD效率高出将近10倍**）。 介绍见：[Bragg Grating full device simulation with EME](https://support.lumerical.com/hc/en-us/articles/360042304334)

但是EME需合理的设置，保证结果的准确。

#### EME 设置：

- 每个光栅周期单元的cell数量要足够，这是保证结果准确的关键前提。（如， period = 317 nm, cell数量需要至少 40 ）
- 边界条件设置为 metal， 比PML速度更快，原因 见：[Starting with metal BCs in FDE simulations](https://support.lumerical.com/hc/en-us/articles/360034901873-Starting-with-metal-BCs-in-FDE-simulations)
- 合理应用对称边界条件提高效率（如若只仿真TE或TM模式，z min可以设置为对称或反对称）
- 计算的EME波长中心越接近所计算的光谱能量中心，结果越准确。（因为EME内置的波长扫描是一个估算。）

## 验证与计算细节 

来自王文康（ wangwenkang@mail.hfut.edu.cn）：

2022-04-05：

> FDTD的Rpeak的近似值为0.4456（估读）
>
> EME仿真100个点的Rpeak近似值为0.4449（估读），存在0.0005左右的差别，增加cell后会更精确，会无限接近FDTD的Rpeak，
>
> 但是cell再增加软件就可能会无响应，个人感觉花费的时间和算力就不值得这么做，因此这次先分为100个cell进行仿真。
>
> 以下是在100个cell，传输模式数量为2的条件下仿真的结果：
>
> 1、没有使用对称的边界条件，
>
> EME使用metal计算边界条件的时间为：6min
>
> EME使用PML计算边界条件的时间为：26min
>
> 由此可见采用metal边界条件仿真时间比PML短
>
> 2、PML条件下：
>
> EME计算出的带宽为34nm
>
> FDTD计算出的带宽为34nm
>
> metal条件下，EME测得的带宽不变
>
> 中心波长一致，均为1544nm-1545nm之间
>
> 3、边界条件全部为PML，未采用对称边界条件时，
>
> FDTD仿真时间：1h
>
> EME仿真时间：26min
>
> 边界条件全部为PML，采用对称边界条件时，
>
> FDTD仿真时间：43min
>
> EME仿真时间：5min
>
> 由此可见采用对称边界条件，EME效果比FDTD更好
>
> 4、中心波长设置为1.61um后，采用EME仿真
>
> 带宽略微变宽一点，大约37nm；
>
> 光谱中心略微左移，依旧为1544-1545nm之间（估读可能存在误差）；
>
> Rpeak变大，为0.493；
>
> 猜测可能是因为波长约长，穿透能力越弱，因此Rpeak会升高，对于较短的波长，穿透能力强，Rpeak会降低。
>
> 之后我将中心波长设置为1.51um进行了验证，
>
> 带宽大约37nm；
>
> 光谱中心略微右移，在1545nm右边；
>
> 总体来说，改变中心波长后光谱中心依旧在1544nm-1545nm之间。
>
> Rpeak变小，为0.41，验证了上边的猜测。



2022-04-04：

> 通过EME这次的仿真结果得到的Rpeak与FDTD一致。
>
> subcell条件、切分cell的数量、传输模式的数量会对实验结果产生影响。
>
> FDTD：Rpeak接近0.45
>
> EME：
>
> 1、更改mesh精度，结果不变
>
> 2、更改subcell条件，由none改为CVCS，反射峰值从略微下降，大概从0.368下降到0.362
>
> 结论：无特殊要求的话就默认subcell条件为none
>
> 3、更改切分cell的数量，可以明显看出反射峰值上升
>
> cell = 10 ；Rpeak = 0.375左右
>
> cell = 20 ；Rpeak = 0.42左右
>
> cell = 30 ；Rpeak = 0.43左右
>
> cell = 40 ；Rpeak = 接近0.45
>
> cell = 50 ，100；Rpeak = 接近0.45，与FDTD仿真结果一致
>
> 结论：cell越多越接近真实值，从数据也可以看出，cell的数量越多，Rpeak越接近FDTD的值。
>
> 4、改变传输模式的数量，可以得：
>
> number of modes for all cell groups = 1，Rpeak = 0.372；
>
> number of modes for all cell groups = 10，Rpeak = 0.372；
>
> number of modes for all cell groups = 20，Rpeak = 0.368；
>
> number of modes for all cell groups = 30，Rpeak = 0.362；
>
> number of modes for all cell groups = 50，Rpeak = 0.35 ；
>
> 结论：随着cell groups的上升，Rpeak逐渐下降，无特殊要求，默认cell groups≤10
>
> 5、改变边界条件，Rpeak不受影响