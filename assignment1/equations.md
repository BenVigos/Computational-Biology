---
geometry:
- margin=1in
---

##  Question 1: The "Living" Ice
### Derivation of $K_{m2}$
$$\begin{align}
v &= \frac{V_{max}S_1S_2}{K_{m1}S_2+K_{m2}S_1+S_1S_2} \\
v(K_{m1}S_2+K_{m2}S_1+S_1S_2) &= V_{max}S_1S_2 \\
K_{m1}S_2+K_{m2}S_1+S_1S_2 &= \frac{V_{max}S_1S_2}{v} \\
K_{m2}S_1 &= \frac{V_{max}S_1S_2}{v} - K_{m1}S_2 - S_1S_2 \\
K_{m2} &= \frac{V_{max}S_2}{v} - \frac{K_{m1}S_2}{S_1} - S_2
\end{align}$$

## Question 2: The Case of the possible Biomass
$$N=\begin{matrix}
CIT \\
ICT \\
AKG \\
SCA \\
SUC \\
FUM \\
MAL \\
OAA \\
X \\


$$Nv=0$$

Gives:
$$\dot{CIT}=2v_1-v_2=0$$
$$\dot{ICT}=v_2-v_3=0$$
$$\dot{AKG}=v_3-v_4=0$$
$$\dot{SCA}=v_4-v_5=0$$
$$\dot{SUC}=v_5-v_6=0$$
$$\dot{FUM}=v_6-v_7=0$$
$$\dot{MAL}=v_7-v_8=0$$
$$\dot{OAA}=v_8-v_1-v_9=0$$
$$\dot{X}=v_9=0$$

From this we get:
$$v_8=v_7=v_6=v_5=v_4=v_3=v_2=2v_1$$

$$\dot{OAA}=v_6-v_1-v_9=v_6-v_1-D=0$$

In order to have biomass conversion to X we must have:
$$\dot{X}=v_9=D> 0$$

We get:
$$v_6-v_1=D$$
Since $D>0$
$$v_6> v_1$$

c.

Irrevirsibility constraint gives $v_1\geq 0$ and MM constraints give $v_6\leq v_{6,max}$.

$$v_6=v_1+D\implies v_1+D\leq v_{6,max}\iff v_1\leq v_{6,max}-D$$

And $v_6>v_1$.

$$v_1\geq0\implies v_1=v_6-D\geq0 \iff v_6 \geq D$$

and $v_6=(V_{max}[Succ])/(K_m+[Succ])-(V_{max}[Fum])/(K_m+[Fum])$.

e. 
