function [data nans]=unnan(data)
% [data nans]=unnan(data) : interpolate at nan, otherwise filtering will kill a lot of data points
% nans contains the positions of the nans in the original data
% eats vectors
    Range=1:size(data,1);
    for i=1:size(data,2);   
      nans=isnan(data(:,i));  

      data(nans,i)=interp1(Range(~nans),data(~nans,i),Range(nans),'pchip');
    end
end