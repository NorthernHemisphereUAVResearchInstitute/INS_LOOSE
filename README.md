# INS_LOOSE
This solution was designed for inertial navigation system of in-pipe inspection. The 15-states extended Kalman filter was used to estimate angle error , speed error,position error,gyroscope drift and accelerometer bias.


*15阶松组合滤波算法*

    惯性导航系统具有自主性、抗干扰性和隐蔽性好，高速连续输出导航参数的优点，但是导航误差随时间累积；卫星导航具有定位和测速精度高，导航结果没有累积误差等优点，但是在弱信号、高动态和强干扰环境下无法跟踪卫星信号，且数据更新率低。将两者组合使用，能实现优势互补。常用的组合方式有松组合、紧组合和深组合。
GNSS/INS组合导航，能够利用GNSS定位结果有效约束INS导航误差，实现优势互补。其中，深组合算法需要自研GNSS接收机，本文不展开研究；紧组合算法鲁棒性低于松组合算法，本代码主要研究松组合算法。

*算法描述*

    系统对陀螺和加速度计数据进行200Hz的采样，以200Hz的速率输出陀螺和加速度计的数据，经过双子样计算，以100Hz的速率输出导航数据。
导航坐标系采用北东地地理坐标系，记为n系，组合系统的捷联坐标系记为b系。采用双字样优化的圆锥效应和划船效应补偿算法，其中现周期均为2，未采用过周期。采用增量进行计算，一个导航计算周期由器件输出周期构成，分别用括号内的1和2表示。计算将相应器件周期内的器件输出扣除零偏，再乘以更新周期，得到增量，再进行误差补偿即可得到速度和姿态的增量，进行导航参数更新。
    陀螺仪测量值经过圆锥误差补偿后，通过旋转矢量法计算姿态更新四元数，利用姿态更新四元数修正方向余弦矩阵，便得到了载体系相对于地理坐标系的转换关系，进一步可以计算得到俯仰角、横滚角和航向角。
利用陀螺仪测量值对加速度计测量值进行划桨效应补偿和旋转效应补偿后，将载体系内的加速度测量值转换到地里坐标系，在地里坐标系内扣除有害加速度以及哥氏力影响后，经过积分可得到位置和速度，完成惯性导航解算。



