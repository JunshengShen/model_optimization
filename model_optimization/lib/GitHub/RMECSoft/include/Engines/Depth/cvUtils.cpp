#include  <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/ml/ml.hpp>
#include <opencv2/highgui/highgui.hpp>

namespace SU{

	/*Draw line to a image with the line normal vector
	*(@lineParameter[0], @lineParameters[1])
	*and a point on line
	*(@lineParameter[2], @lineParameters[3])
	*equation:  ax + bx = c
	*           a = nx = @lineParameter[0]
	*           b = ny = @lineParameter[1]
	*/
	void drawLine(cv::Mat mat, int x1, int x2,std::vector<double> lineParameters, int colorIdx, int lineWidth)
	{
		double c = lineParameters[0] * lineParameters[2] + lineParameters[1] * lineParameters[3];
		double k = lineParameters[0]/lineParameters[1];
		c = c/lineParameters[1];

		int y1 = (int)(c - k * x1);
		int y2 = (int)(c - k * x2);

		if (colorIdx == -1)
		{
			cv::line(mat, cv::Point(x1, y1), cv::Point(x2,y2), cv::Scalar(0,255,0), lineWidth);
			return;
		}
		cv::line(mat, cv::Point(x1, y1), cv::Point(x2,y2), cv::Scalar(0,0,255), lineWidth);

	}

	/*Draw line on median-filtered buffer
	*@lenFilter   length of median filter buffer
	*(@lineParameter[0], @lineParameters[1])
	*and a point on line
	*(@lineParameter[2], @lineParameters[3])
	*equation:  ax + by = c
	*           a = nx = @lineParameter[0]
	*           b = ny = @lineParameter[1]
	*/
	void drawLineScale(cv::Mat mat, int lenFilter,int x1, int x2,std::vector<double> lineParameters, int colorIdx, int lineWidth)
	{
		if (lineParameters.empty())
		{
			return;
		}
		double c = lineParameters[0] * lineParameters[2] + lineParameters[1] * lineParameters[3];
		double k = lineParameters[0]/lineParameters[1];		
		c = c/lineParameters[1];

		int y1 = (int)(c - k * x1);
		int y2 = (int)(c - k * x2);

		if (colorIdx == -1)
		{
			cv::line(mat, cv::Point(x1/lenFilter, y1), cv::Point(x2/lenFilter,y2), cv::Scalar(0,255,0), lineWidth);
			return;
		}
		cv::line(mat, cv::Point(x1/lenFilter, y1), cv::Point(x2/lenFilter,y2), cv::Scalar(0,0,255), lineWidth);

	}

}
