#version 3.6;
global_settings{assumed_gamma 1.0}
#include "colors.inc"
#include "textures.inc"
// camera -----------------------------------------------------------
camera {location <50, 50, 100>
look_at  <50, 50,  0> }
// sun ---------------------------------------------------------------
light_source{<0,0,2500> color White}
// sky -------------------------------------------------------------
//background
background{White}

//========== grid=======================
union{
#fopen grid_file "GRID.dat" read
#while(defined(grid_file))
    #read (grid_file, t1, t2, t3, t4)
    cylinder{<t1, t2, 0>,
        <t3, t4, 0>,
        0.1 //grid w
        pigment{Black transmit 0.9}
        finish{ambient.1}
    }
#end
#fclose grid_file
}

//molecules
union{
#fopen atom_file "POV.DATA" read
#while(defined(atom_file))
    #read (atom_file, c1,c2)
        polygon{4,
        <c1-0.5, c2-0.5>, <c1-0.5,c2+2+0.5 >, <c1+2+0.5, c2+2+0.5>, <c1+2+0.5, c2-0.5>
        texture { pigment { Red }
        }
    }
#end
#fclose atom_file
}


//f-h
union{
#fopen deH_file "fH.DATA" read
#while(defined(deH_file))
    #read (deH_file, c1,c2)
        sphere{<c1, c2, 0>,0.6
            pigment {Yellow transmit 0.2}
            finish{ambient.3}
            }
#end
#fclose deH_file
}


//bond
union{
#fopen bond_file "bond.DATA" read
#while(defined(bond_file))
#read (bond_file,c1,c2,c3,c4)
    cylinder{<c1, c2, 0>, <c3, c4, 0>,0.4
        pigment {Blue}
        finish{ambient.3}
    }
#end
#fclose bond_file
}

//Hydrogen
union{
#fopen hydrogen_file "HYDROGEN.DATA" read
#while(defined(hydrogen_file))
#read (hydrogen_file,t1,t2)
sphere{<t1, t2, 0>,0.6
    pigment {Yellow transmit 0.2}
    finish{ambient.3}
    }
#end
#fclose hydrogen_file
}
