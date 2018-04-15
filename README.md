# RayTracingPPM
A simple ray tracing project which generates PPM file.

——————————————————————

更新日志：2017-07-12

加入了光源、光圈模糊、和并行处理。
1、关于光源，效果参见Example7，或者选择场景的时候输入1.
2、关于光圈模糊，效果参见Example8。或在程序进行相关询问的时候输入1.
3、关于并行处理……我使用G++编译的，编译的时候加了 -fopenmp……我也不知道VC++怎么加这个
4、Example5，和Example6是之前图片的高配版本

——————————————————————

1、关于编译，因为就一个.cpp和一个.h 随便找个编译器就能编译，需要支持 C++11 ，为了加速运行最好加-O3。
2、所有代码都是手敲的……没有复制粘贴，很多抄了书但也是手打的……有一部分数学代码是上学期写的。
3、因为我电脑太差了！！！所以设了分辨率上限是1000x1000，代码里可以改……抗锯齿的采样次数也写死成75了， 代码里有一个叫ns的变量可以改…… 
4、目前看来循环次数为75，分辨率1000x1000，十分钟之内应该可以跑出来。
5、Example123分别是同一个场景在不同参数下的表现 ，Example4是一个随机场景。参数都不高……高了跑不动…… 
6、让输入场景的时候，如果你输入0，会出现样例场景，如果输入其他，会出现随机场景。
7、代码很差没有鲁棒性……如果输数字的时候输入了字符……Emmmm我也不知道会咋样……
8、太菜了……阴影光源对焦模糊都没加……  
啊 如果不知道输入什么，
第一次输入50
第二次输入800
第三次输入800
第四次输入0或者1
就行了…… 
