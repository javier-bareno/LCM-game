#include "colors.inc"           
#include "textures.inc"
 
 background { color White }    
 
 
 
  /**********************************************
 ** Camera and lights                              
 **********************************************/         
 
#declare camdist = 200;
                                                   
camera {   
        perspective orthographic
        location <0, 0, camdist>
        look_at <0, 0, 0>
}    


light_source { <0, 0, 10000> color White}
light_source { <-5000, 00, 10000> color White}   
//light_source { <000, 00, 1.5*camdist> color White}                                                  
//light_source { <-600, 00, -camdist> color White}
//light_source { <600, 00, -camdist> color White}


 /**********************************************
 ** Macros for drawing atoms
 **********************************************/  
 
#declare RTiN = 1.35;       
#declare Li_col = OldGold;
#declare Co_col =  MediumBlue;
#declare Mn_col = Sienna;

plane { <0,0,1>, -2 pigment {color 0.20* Yellow+ 0.8* White}}

 #macro PutAtom_Li(X,Y,Z)
/* sphere{<X,Y,Z>,RTiN 
 texture{
   pigment{color Li_col     }
          // transmit 0.65} 
   finish { phong 0.7 }   
  }}                  */
#end   
                

#macro PutAtom_Mn(X,Y,Z)
 sphere{<X,Y,Z>,RTiN 
 
 texture { 
   pigment{color Mn_col} 
   finish { phong 0.7 }
  }} 
#end
     
     
#macro PutAtom_Co(X,Y,Z)
 sphere{<X,Y,Z>,RTiN 
 texture { 
   pigment{color Co_col} 
   finish { phong 0.7 }
  }}  
#end 
 



union {

        #include "..\pickle_3hex.pov"   
        
        
        // ad rotate and anim code here  
        
        //translate <-5000, -3600, 0>
       //rotate <0, 0, 90> 
        
        
        #local ang = clock * 90;
        rotate <-ang, 0, 0> 
        
}
