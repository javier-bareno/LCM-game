#include "colors.inc"           
#include "textures.inc"
 
 background { color White }    
 
 
 
  /**********************************************
 ** Camera and lights                              
 **********************************************/         
 
#declare camdist = 150;
                                                   
camera {   
        perspective orthographic
        location <0, 0, camdist>
        look_at <0, 0, 0>
}    


light_source { <0, 0, 100> color White}
light_source { <-50, 00, 100> color White}   
//light_source { <000, 00, 1.5*camdist> color White}                                                  
//light_source { <-600, 00, -camdist> color White}
//light_source { <600, 00, -camdist> color White}


 /**********************************************
 ** Macros for drawing atoms
 **********************************************/  
 
#declare RTiN = 1.35;       
#declare Li_col = OldGold;
#declare Co_col = MediumBlue;
#declare Mn_col = SteelBlue;

plane { <0,0,1>, -2 pigment {color Bronze}}
  
                
 #macro PutAtom_(a,b,X,Y,Z)
 sphere{<X,Y,Z>,RTiN    
 texture{
   pigment{#switch(a)
             #case(0)
               #declare s = b/1300;
               color s * Blue + (1-s) * White  }
             #break
             #case(1)  
               #declare s = b/1300;
               color s * Green + (1-s) * White  }  
             #break
           #end
           //transmit 0.65} 
   finish { phong 0.7 } 
  }}  
#end 

 



union {

        #include "test_4.pov"   
        
        
        // ad rotate and anim code here  
        
       translate -48.5*<0.5,0.866,0> 
       rotate <0,180, 0> 
        
        
        
        //rotate <0, 0, 90> 
        
}
