# Machine Learning Code
The code contained in this folder is for the ML chapter of work within the thesis. Several files are included.

Image Rekognition should have the files uploaded into an S3 bucket, then have the code from Sagemaker process the data into the slices, 
and using the Image Rekognition interface, a model should be trained upon the data. Going back to the sagemaker model, it is possible 
to deploy the Rekognition model by copying across the relevant information about the trained model's Amazon Resource Name (ARN).

## WARNING: When done using the model, ensure the model is stopped on the Rekognition interface to ensure costs are not incurred when not in use.
