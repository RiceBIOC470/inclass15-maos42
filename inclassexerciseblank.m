%% step 1: write a few lines of code or use FIJI to separately save the
% nuclear channel of the image Colony1.tif for segmentation in Ilastik
%ih10@rice.edu

img=('48hColony1.tif');
ishow=imread(img);
reader=bfGetReader(img);

z=reader.getSizeZ;
time=reader.getSizeT;
channel=reader.getSizeC;

iplane=reader.getIndex(z-1,channel-4,time-1)+1;
DAPI=bfGetPlane(reader,iplane);
imshow(DAPI,[]);

imwrite(DAPI,'48hColony1channel1withMATLAB.tif');


%% step 2: train a classifier on the nuclei
% try to get the get nuclei completely but separe them where you can
% save as both simple segmentation and probabilities

%% step 3: use h5read to read your Ilastik simple segmentation
% and display the binary masks produced by Ilastik 

% (datasetname = '/exported_data')
% Ilastik has the image transposed relative to matlab
% values are integers corresponding to segmentation classes you defined,
% figure out which value corresponds to nuclei

segs=h5read('48hColony1_DAPI_Simple Segmentation.h5', '/exported_data');
probs=h5read('48hColony1_DAPI_Probabilities.h5', '/exported_data');

h5disp('48hColony1_DAPI_Simple Segmentation.h5', '/exported_data');
h5disp('48hColony1_DAPI_Probabilities.h5', '/exported_data');

imshow(squeeze(segs==2)) %squeeze gets rid of useless column. the value 2, 
%its because on ilastik you defined channel 2 as the nucleus green
img2=squeeze(segs==2);

%% step 3.1: show segmentation as overlay on raw data

imshow(ishow,[]); hold on; imshow(img2);
imshowpair(ishow, img2, 'montage');

img3=imfuse(ishow,img2);
imshow(img3);
%% step 4: visualize the connected components using label2rgb
% probably a lot of nuclei will be connected into large objects

img4=label2rgb(img2); %on the segmentation?
imshow(img4,[]);

%% step 5: use h5read to read your Ilastik probabilities and visualize
% it will have a channel for each segmentation class you defined

newprobs=zeros(2301,2);
%newprobs(:,1:2)=probs(:,:,2:3); %doesn't work :( 

imgprob=squeeze(probs(2,:,:));
figure (3); imshow(imgprob);

%% step 6: threshold probabilities to separate nuclei better

averageprob=mean(mean(imgprob));
stdevprob=std(std(imgprob));

newprob=imgprob>averageprob+2*stdevprob;
figure (4); imshow(newprob);

%% step 7: watershed to fill in the original segmentation (~hysteresis threshold)

%this one works
out=~imdilate(img2,strel('disk',3));
basin=imcomplement(bwdist(out));
basin=imimposemin(basin,newprob|out);
L=watershed(basin);
L(~img2)=0;
hope=label2rgb(L,'jet','k','shuffle');
figure (5); imshow(hope);


%% step 8: perform hysteresis thresholding in Ilastik and compare the results
% explain the differences

hyst=h5read('hysteresis.h5', '/exported_data');
imghyst=squeeze(hyst);
RGB=label2rgb(imghyst,'jet','k','shuffle');
figure (6); imshow(RGB,[]);

%Miguel Angel: The resulting hysteresis from ilastik has more cells
%separate, however it does contain greater empty areas, and the nuclei have holes in them. 
%compared to the watershed in part 7 even though it seems all cells are
%conserved, they are not entirely separte, as we can see the big areas from
%colors,meaning a better separtion is needed. 



%% step 9: clean up the results more if you have time 
% using bwmorph, imopen, imclose etc

BW2=bwmorph(img2,'open');
BW2=bwmorph(img2,'thin');
imshow(BW2);

out2=~imdilate(BW2,strel('disk',3));
basin2=imcomplement(bwdist(out2));
basin2=imimposemin(basin2,newprob|out2);
L2=watershed(basin2);
L2(~BW2)=0;
hope2=label2rgb(L2,'jet','k','shuffle');
figure (7); imshow(hope2);
