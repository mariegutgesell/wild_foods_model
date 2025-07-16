= Functional Response Model for Omnivory and Habitat Preference

== Habitat Preference and Foraging Scale

Habitat preference (Wi) is density-dependent, with fixed preference `w`:
#math Wi = \frac{w C_1}{w C_1 + (1 - w) C_2}

Foraging scale (S) depends on this preference and a scaling parameter `Q`:
#math S = Q + (1 - Q) Wi

When `Q = 0`, foraging is exclusive (single patch). When `Q = 1`, consumers perceive prey as well-mixed.

== Degree of Omnivory

Defined as density-dependent preference for basal vs. intermediate prey (McCann et al. 2005):

Patch 1:
#math o_1 = \frac{o e a R_1}{o e a R_1 + (1 - o) e a C_1}

Patch 2:
#math o_2 = \frac{o e a R_2}{o e a R_2 + (1 - o) e a C_2}

== Active Omnivory Preference (Gutgesell et al. 2022)

Patch 1:
#math o_1^\text{pref} = \frac{o R_1}{o R_1 + (1 - o) C_1}

Patch 2:
#math o_2^\text{pref} = \frac{o R_2}{o R_2 + (1 - o) C_2}


== Grocery Preference

Supplementary foraging on groceries `G` depends on a switching parameter `H`:
#math H_G = \frac{H G}{H G + (1 - H)(R_1 + R_2 + C_1 + C_2)}

== Basic Functional Response (Type II)



== Predator Functional Responses



== Consumer Functional Response



== Differential Equations

Let state vector be:
#math u = (R_1, R_2, C_1, C_2, P)

Then the ODEs are:
#math \dot{R_1} = r R_1 \left(1 - \frac{R_1}{K}\right) - C_1 f_{R,C1} - P f_{R1,P}
#math \dot{R_2} = r R_2 \left(1 - \frac{R_2}{K}\right) - C_2 f_{R,C2} - P f_{R2,P}
#math \dot{C_1} = e C_1 f_{R,C1} - P f_{C1,P} - m C_1
#math \dot{C_2} = e C_2 f_{R,C2} - P f_{C2,P} - m C_2
#math \dot{P} = 0 \quad \text{(fixed for equilibrium analysis)}


#math Wi = frac(w * C.sub(1), w * C.sub(1) + (1 - w) * C.sub(2))
#math S = Q + (1 - Q) * Wi
#math o.sub(1) = frac(o * e * a * R.sub(1), o * e * a * R.sub(1) + (1 - o) * e * a * C.sub(1))
#math o.sub(2) = frac(o * e * a * R.sub(2), o * e * a * R.sub(2) + (1 - o) * e * a * C.sub(2))
#math o.sub(1).sup(pref) = frac(o * R.sub(1), o * R.sub(1) + (1 - o) * C.sub(1))
#math o.sub(2).sup(pref) = frac(o * R.sub(2), o * R.sub(2) + (1 - o) * C.sub(2))
#math o.sup(pref) = o.sub(fixed)
#math H_G = frac(H * G, H * G + (1 - H) * (R.sub(1) + R.sub(2) + C.sub(1) + C.sub(2)))
#math f(R) = frac(a * R, 1 + a * h * R)
#math f.sub(R1,P) = frac((1 - H) * S * a * o.sub(1) * R.sub(1), 1 + S * a * h * o.sub(1) * R.sub(1) + (1 - S) * a * h * o.sub(2) * R.sub(2) + S * a * h * (1 - o.sub(1)) * C.sub(1) + (1 - S) * a * h * (1 - o.sub(2)) * C.sub(2) + H * a * h * G)
#math f.sub(R2,P) = "same numerator with (1 - S) * a * o_2 * R_2, same denominator as above"
#math f.sub(C1,P) = "similar structure, replace o with (1 - o)"
#math f.sub(C2,P) = "similar structure, with (1 - S) and C_2"
#math f.sub(G,P) = frac(H * a * G, same denominator as above)
#math f.sub(G,P) = G * H
#math f.sub(R,C1) = frac(a * R.sub(1), 1 + a * h * R.sub(1))
#math f.sub(R,C2) = frac(a * R.sub(2), 1 + a * h * R.sub(2))
#math dot(R.sub(1)) = r * R.sub(1) * (1 - R.sub(1)/K) - C.sub(1) * f.sub(R,C1) - P * f.sub(R1,P)
#math dot(R.sub(2)) = r * R.sub(2) * (1 - R.sub(2)/K) - C.sub(2) * f.sub(R,C2) - P * f.sub(R2,P)
#math dot(C.sub(1)) = e * C.sub(1) * f.sub(R,C1) - P * f.sub(C1,P) - m * C.sub(1)
#math dot(C.sub(2)) = e * C.sub(2) * f.sub(R,C2) - P * f.sub(C2,P) - m * C.sub(2)
#math dot(P) = 0

