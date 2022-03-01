# Data set

## MMR-MK instance format
> "Furini, F., Iori, M., Martello, S., & Yagiura, M. (2015). Heuristic and exact algorithms for the interval min–max regret knapsack problem. INFORMS Journal on Computing, 27(2), 392-405."

There are 486 data files.
  
An MMR-KP instance denoted by **_v-ww-xx-yy-zz_** indicates that it is generated based on type **_v_** KP instance with **_ww_** items, **_xx = 100 g_**, and the max **_zz_**%-uncertainty on profit intervals.

The format for each of these data files is:
~~~
number of items (n)
resource capacity (b)
resource consumed in selecting item j (j=1,...,n)
minimum profit of item j (j=1,...,n)
maximum profit of item j (j=1,...,n)
~~~

## MMR-MKP instance format
> "Wu, W., Iori, M., Martello, S., & Yagiura, M. (2020). An Iterated Dual Substitution Approach for Binary Integer Programming Problems under the Min-Max Regret Criterion. arXiv preprint arXiv:2012.07530."

There are 270 data files.

An MMR-MKP instance denoted by **_wwxxxyy-zz_** is generated from the **_zz_**-th MKP instance with **_ww_**-dimensions, **_xxx_**-items and the max **_yy_**%-uncertainty on profit intervals.

The format for each of these data files is:
~~~
number of dimensions (m)	number of items (n)
minimum profit of item j (j=1,...,n)
maximum profit of item j (j=1,...,n)
for each dimension i (i=1,...,m) in turn:
    resource consumed in selecting item j in dimension i (j=1,...,n)
resource capacity of dimension i (i=1,...,m)
~~~

## MMR-SCP instance format
> "Pereira, J., & Averbakh, I. (2013). The robust set covering problem with interval data. Annals of Operations Research, 207(1), 217-235."

There are 225 data files.

An MMR-SCP instance denoted by **_Bxyyzz_** indicates a Type-**_B_** instance whose corresponding SCP instance is the **_yy_**-th instance in family **_x_** from the OR-Library with the max **_zz_**%-uncertainty on cost intervals,
while **_Mxyy-z_** (or **_Kxyy-z_**) stands for the **_z_**_th Type-**_M_** (or Type-**_K_**) instance whose corresponding SCP instance is the **_yy_**-th instance in family **_x_**.

The format for each of these data files is:
~~~
number of agents (m)
number of jobs (n)
for each agent i (i=1,...,m) in turn:
    minimum cost of allocating job j to agent i (j=1,...,n)
for each agent i (i=1,...,m) in turn:
    maximum cost of allocating job j to agent i (j=1,...,n)
for each agent i (i=1,...,m) in turn:
    resource consumed in allocating job j to agent i (j=1,...,n)
resource capacity of agent i (i=1,...,m)
~~~

## MMR-GAP instance format
> "Wu, W., Iori, M., Martello, S., & Yagiura, M. (2018). Exact and heuristic algorithms for the interval min-max regret generalized assignment problem. Computers & Industrial Engineering, 125, 98-110."

There are 300 data files and 60 files for each type.
  
An MMR-GAP instance denoted by **_Txxyyzz-i_** indicates the **_i_**-th instance in a group of problems of type **_T_** with **_xx_**-agents, **_yy_**-jobs and the max **_zz_**%-uncertainty on cost intervals.

The format for each of these data files is:
~~~
number of agents (m)
number of jobs (n)
for each agent i (i=1,...,m) in turn:
    minimum cost of allocating job j to agent i (j=1,...,n)
for each agent i (i=1,...,m) in turn:
    maximum cost of allocating job j to agent i (j=1,...,n)
for each agent i (i=1,...,m) in turn:
    resource consumed in allocating job j to agent i (j=1,...,n)
resource capacity of agent i (i=1,...,m)
~~~