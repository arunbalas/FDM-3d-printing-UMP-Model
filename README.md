## A System and Architecture of Fused Deposition Modeling - Unit Manufacturing Process (FDM-UMP)

In today’s competitive world economy, the manufacturing and design engineers face the challenge of manufacturing components rapidly to meet customer requirements and achieve competitive edge. Additive manufacturing provides an efficient method to build complex products or prototypes to minimize the design and cycle time. Fused Deposition Modeling (FDM) is an additive manufacturing process used to build prototypes using variety of materials. The build and support materials are extruded as a semi-molten filament through the extrusion head and deposited layer by layer to construct prototypes directly from 3D CAD model. This technology is increasingly used for customized products, conceptual models and finds its applications in many fields of engineering and industry like aerospace, automotive products, dentistry and medical implants etc.
<p align="center">
<img src="https://github.com/arunbalas/FDM-UMP/blob/master/IMG_2093.JPG" width="700" height="500">
</p>
When compared to traditional manufacturing processes like milling, drilling, etc., there are not many mathematical abstractions available for characterizing the FDM process. Majority of the literature is focused upon developing models for specific purposes with limited number of parameters, giving insight on how the process behaves with respect to changing parameters and methodology, thereby finding the optimum levels of parameters. In these cases, only the methodology is considered to be useful, as it hinders the use of developed model for different kinds of machines operating under the same FDM principle. So, developing a general abstraction of FDM process is of high importance to the growing user community which helps in saving cost, time and environment, with a major challenge being the limited amount of data and resources. So, in this paper, FDM-UMP model is developed and validated by fabricating a simple component shown in below figure with the aid of Stratasys Uprint SE Plus printer. The video demo of FDM-UMP model using python is available at this link: https://youtu.be/IaSZnU1Fvy8

<p align="center">
<img src= "https://github.com/arunbalas/FDM-UMP/blob/master/Turbine%20processed%20data.JPG">
Figure: A Sample CAD model used for process parameter estimation.
</p>

<p align="center">
	
![ump](https://github.com/arunbalas/FDM-UMP/blob/master/Graphical%20UMP.jpg)
Figure: Graphical UMP model for FDM (3D Printing) process.
</p>


### REPORT Folder:
This folder consists of the PDF document [FDM-UMP_Arun.pdf] containing all the necessary information of FDM UMP.
	


### CODE Folder:
This folder consists of input.json and Python.py file. You can modify the input process parameters in the json file and run the python code to get the output process parameters (also in json format).


### REFERENCES Folder:
This folder consists of all the main references used for this project.
