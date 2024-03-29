# RKF45-Runge_Kutta_3_order

the program solves a system of differential equations

$\displaystyle\frac{\mathrm d x_1}{\mathrm d t} = -14x_1 + 13x_2 + cos(1 + t); \quad \frac{\mathrm d x_2}{\mathrm d t} = 20x_1 - 30x_2 + arctg(1 + t^2)$

$x_1(0) = 2, \quad x_2(0) = 0.5, \quad t \in [0; 1.5]$;  

in $h_{print} = 0.075$ according to the RKF45 program with $EPS = 0.001$ and the order 3 Runge-Kutta method

$z_{n+1} = z_n + \displaystyle\frac{2k_1 + 3k_2 + 4k_3}{9}, \quad k_1 = hf(t_n, z_n), \quad k_2 = hf(t_n + \displaystyle\frac{h}{2}, z_n + \displaystyle\frac{k_1}{2}), \quad k_3 = hf(t_n + \displaystyle\frac{3h}{4}, z_n + \displaystyle\frac{3k_2}{4})$

in $h_{int} = 0.075$ and $h_{int} = 0.00075$

The first step gives a qualitatively incorrect solution, since it does not satisfy the stability condition

_`Forsythe.h` is a librarian program_
