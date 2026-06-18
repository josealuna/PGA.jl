module PGA

using Grassmann


see(a,b) = a + b
export see

export IIa, IIb
#@basis D"1,1,1,0" # PGA 3D con métrica D"1,1,1,0"

"""
    IIa(P, Ixx, Iyy, Izz, m)

Calcula el bivector de velocidad (twist) `B` a partir del bivector de momento (wrench) `P` 
en PGA 3D (métrica D"1,1,1,0"), aplicando el inverso del operador de inercia y el mapa dual.

# Argumentos
- `P`: Multivector de Grassmann que representa el momento.
- `Ixx`, `Iyy`, `Izz`: Momentos principales de inercia (escalares).
- `m`: Masa total del sólido rígido (escalar).

# Notas
Utiliza el acoplamiento estructural del operador `!` (complemento) verificado para la base 
ordenada (v1, v2, v3, v4=e0), corrigiendo las inversiones de signo producto de las permutaciones.
"""
function IIa(P,Ixx, Iyy, Izz, m)

    # Simple function to compute the inertia of a rigid body in PGA 3D given the momentum and mass properties
    compOf(base,inertia) = (P[basis(!base)] / inertia)*base
    # Return the inertia as a linear combination of the basis elements
    compOf(v12,Izz) +
    -compOf(v13,Iyy) + 
    compOf(v23,Ixx) + 
    compOf(v14,m)   +
    -compOf(v24,m)   +
    compOf(v34,m)

end


"""
    IIb(P, Ixx, Iyy, Izz, m)

Calcula el bivector de velocidad (twist) `B` a partir del bivector de momento (wrench) `P` 
en PGA 3D (métrica D"1,1,1,0"), aplicando el inverso del operador de inercia y el mapa dual.

# Argumentos
- `P`: Multivector de Grassmann que representa el momento.
- `Ixx`, `Iyy`, `Izz`: Momentos principales de inercia (escalares).
- `m`: Masa total del sólido rígido (escalar).

# Notas
Utiliza el acoplamiento estructural del operador `!` (complemento) verificado para la base 
ordenada (v1, v2, v3, v4=e0), corrigiendo las inversiones de signo producto de las permutaciones.
"""
function IIb(P, Ixx, Iyy, Izz, m)
    # Mapeo directo corrigiendo los signos del dual combinatorio '!' de tu laptop:
    # v12 <-> v34 (sin cambio de signo)
    # v13 <-> v24 (cambia de signo en el operador !)
    # v23 <-> v14 (sin cambio de signo)
    
    B = (P[v34] / Izz) * v12 +
        (-P[v24] / Iyy) * v13 +  # Ajuste de signo por el dual de v13 -> -v24
        (P[v14] / Ixx) * v23 +
        (P[v23] / m)   * v14 +
        (-P[v13] / m)  * v24 +  # Ajuste de signo por el dual de v24 -> -v13
        (P[v12] / m)   * v34
        
    return B
end


end