# INS_LOOSE
This solution was designed for inertial navigation system of in-pipe inspection. The 15-states extended Kalman filter was used to estimate angle error , speed error,position error,gyroscope drift and accelerometer bias.


惯性导航系统具有自主性、抗干扰性和隐蔽性好，高速连续输出导航参数的优点，但是导航误差随时间累积；卫星导航具有定位和测速精度高，导航结果没有累积误差等优点，但是在弱信号、高动态和强干扰环境下无法跟踪卫星信号，且数据更新率低。将两者组合使用，能实现优势互补。常用的组合方式有松组合、紧组合和深组合。
GNSS/INS组合导航，能够利用GNSS定位结果有效约束INS导航误差，实现优势互补。其中，深组合算法需要自研GNSS接收机，本文不展开研究；紧组合算法鲁棒性低于松组合算法，本代码主要研究松组合算法
