/*---------------------------------------------------------------------
*
* Copyright © 2016  Minsi Chen
* E-mail: m.chen@derby.ac.uk
*
* The source is written for the Graphics I and II modules. You are free
* to use and extend the functionality. The code provided here is functional
* however the author does not guarantee its performance.
---------------------------------------------------------------------*/
#include <algorithm>
#include <math.h>

#include "Rasterizer.h"

Rasterizer::Rasterizer(void)
{
	mFramebuffer = NULL;
	mScanlineLUT = NULL;
}

void Rasterizer::ClearScanlineLUT()
{
	Scanline *pScanline = mScanlineLUT;

	for (int y = 0; y < mHeight; y++)
	{
		(pScanline + y)->clear();
		(pScanline + y)->shrink_to_fit();
	}
}

unsigned int Rasterizer::ComputeOutCode(const Vector2 & p, const ClipRect& clipRect)
{
	unsigned int CENTRE = 0x0;
	unsigned int LEFT = 0x1;
	unsigned int RIGHT = 0x1 << 1;
	unsigned int BOTTOM = 0x1 << 2;
	unsigned int TOP = 0x1 << 3;
	unsigned int outcode = CENTRE;
	
	if (p[0] < clipRect.left)
		outcode |= LEFT;
	else if (p[0] >= clipRect.right)
		outcode |= RIGHT;

	if (p[1] < clipRect.bottom)
		outcode |= BOTTOM;
	else if (p[1] >= clipRect.top)
		outcode |= TOP;

	return outcode;
}

bool Rasterizer::ClipLine(const Vertex2d & v1, const Vertex2d & v2, const ClipRect& clipRect, Vector2 & outP1, Vector2 & outP2)
{
	//TODO: EXTRA This is not directly prescribed as an assignment exercise. 
	//However, if you want to create an efficient and robust rasteriser, clipping is a usefull addition.
	//The following code is the starting point of the Cohen-Sutherland clipping algorithm.
	//If you complete its implementation, you can test it by calling prior to calling any DrawLine2D .

	const Vector2 p1 = v1.position;
	const Vector2 p2 = v2.position;
	unsigned int outcode1 = ComputeOutCode(p1, clipRect);
	unsigned int outcode2 = ComputeOutCode(p2, clipRect);

	outP1 = p1;
	outP2 = p2;
	
	bool draw = false;

	return true;
}

void Rasterizer::WriteRGBAToFramebuffer(int x, int y, const Colour4 & colour)
{
	PixelRGBA *pixel = mFramebuffer->GetBuffer();
	Colour4 newColour = colour;
	if (x > 0 && x < mWidth && y > 0 && y < mHeight - 1)
	{
		if (mBlendMode == ALPHA_BLEND)
		{
			Colour4 old = pixel[y*mWidth + x];
			Colour4 temp(colour * colour[3] + old * (1.0 - colour[3]));
			newColour.SetVector(temp[0], temp[1], temp[2], temp[3]);
		}

		pixel[y*mWidth + x] = newColour;
	}
}

Rasterizer::Rasterizer(int width, int height)
{
	//Initialise the rasterizer to its initial state
	mFramebuffer = new Framebuffer(width, height);
	mScanlineLUT = new Scanline[height];
	mWidth = width;
	mHeight = height;

	mBGColour.SetVector(0.0, 0.0, 0.0, 1.0);	//default bg colour is black
	mFGColour.SetVector(1.0, 1.0, 1.0, 1.0);    //default fg colour is white

	mGeometryMode = LINE;
	mFillMode = UNFILLED;
	mBlendMode = NO_BLEND;

	SetClipRectangle(0, mWidth, 0, mHeight);
}

Rasterizer::~Rasterizer()
{
	delete mFramebuffer;
	delete[] mScanlineLUT;
}

void Rasterizer::Clear(const Colour4& colour)
{
	PixelRGBA *pixel = mFramebuffer->GetBuffer();

	SetBGColour(colour);

	int size = mWidth*mHeight;
	
	for(int i = 0; i < size; i++)
	{
		//fill all pixels in the framebuffer with background colour
		*(pixel + i) = mBGColour;
	}
}

void Rasterizer::DrawPoint2D(const Vector2& pt, int size)
{
	int x = pt[0];
	int y = pt[1];
	
	WriteRGBAToFramebuffer(x, y, mFGColour);
}

void Rasterizer::DrawLine2D(const Vertex2d & v1, const Vertex2d & v2, int thickness)
{
	//The following code is basic Bresenham's line drawing algorithm.
	//The current implementation is only capable of rasterise a line in the first octant, where dy < dx and dy/dx >= 0
	//See if you want to read ahead https://www.cs.helsinki.fi/group/goa/mallinnus/lines/bresenh.html

	//TODO:
	//Ex 1.1 Complete the implementation of Rasterizer::DrawLine2D method. 
	//This method currently consists of a partially implemented Bresenham algorithm.
	//You must extend its implementation so that it is capable of drawing 2D lines with arbitrary gradient(slope).
	//Use Test 1 (Press F1) to test your implementation

	//Ex 1.2 Extend the implementation of Rasterizer::DrawLine2D so that it is capable of drawing lines based on a given thickness.
	//Note: The thickness is passed as an argument int thickness.
	//Use Test 2 (Press F2) to test your implementation

	//Ex 1.3 Extend the implementation of Rasterizer::DrawLine2D so that it is capable of interpolating colour across a line when each end-point has different colours.
	//Note: The variable mFillMode indicates if the fill mode is set to INTERPOLATED_FILL. 
	//The colour of each point should be linearly interpolated using the colours of v1 and v2.
	//Use Test 2 (Press F2) to test your implementation

	Vector2 pt1 = v1.position;
	Vector2 pt2 = v2.position;


	bool swap_vertices = pt1[0] > pt2[0];

	if (swap_vertices)
	{
		pt1 = v2.position;
		pt2 = v1.position;
	}

	int dx = pt2[0] - pt1[0];
	int dy = pt2[1] - pt1[1];
	
	int epsilon = 0;

	bool negative_slope = dy < 0;
	bool swap_xy = (abs(dx) < abs(dy));

	int x = pt1[0];
	int y = pt1[1];
	if (negative_slope && !swap_xy)
	{
		y = -y;
		dy = -dy;
	}
	int ex = pt2[0];

	//If gradient is too steep
	if (swap_xy)
	{
		//Octant 7, reflect around x
		if (negative_slope && !swap_vertices)
		{
			y = -y;
			dy = -dy;
		}

		int temp(x);
		x = y;
		y = temp;

		int tempdx(dx);
		dx = dy;
		dy = tempdx;

		ex = pt2[1];

		//Octant 3, reflect around x
		if (negative_slope && swap_vertices)
		{
			x = -x;
			dx = -dx;
		}
	}

	while (swap_xy && negative_slope ? ex <= (x < 0 ? abs(x) : -x) : x <= ex)
	{
		Colour4 colour = v1.colour;
		
		Vector2 temp;
		if (swap_xy)
		{
			temp[0] = y;
			temp[1] = negative_slope ? -x : x;
		}
		else
		{
			temp[0] = x;
			temp[1] = negative_slope ? -y : y;
		}

		if (mFillMode == INTERPOLATED_FILLED)
		{
			float t = abs(pt1.Norm() - temp.Norm()) / abs(pt2.Norm() - pt1.Norm());
			colour = (swap_vertices ? v2.colour : v1.colour) * t + (swap_vertices ? v1.colour : v2.colour) * (1 - t);
		}
		SetFGColour(colour);

		temp[1] += thickness % 2 == 0 ? thickness / 2 - 1 : (thickness - 1) / 2;

		for (int i = 0; i < ((thickness % 2 == 0) ? (thickness / 2) - 1 : (thickness - 1) / 2); i++)
		{
			temp[1] -= 1;
			DrawPoint2D(temp);
		}
		DrawPoint2D(temp); 

		for (int i = 0; i < ((thickness % 2 == 0) ? (thickness / 2) : (thickness - 1) / 2); i++)
		{
			temp[1] += 1;
			DrawPoint2D(temp);
		}
		
		epsilon += dy;
		
		if ((epsilon << 1) >= dx)
		{
			y++;
			epsilon -= dx;
		}

		x++;
	}
}

void Rasterizer::DrawUnfilledPolygon2D(const Vertex2d * vertices, int count)
{
	//TODO:
	//Ex 2.1 Implement the Rasterizer::DrawUnfilledPolygon2D method so that it is capable of drawing an unfilled polygon, i.e. only the Edges of a polygon are rasterised. 
	//Please note, in order to complete this exercise, you must first complete Ex1.1 since DrawLine2D method is reusable here.
	//Note: The Edges of a given polygon can be found by conntecting two adjacent vertices in the vertices array.
	//Use Test 3 (Press F3) to test your solution.

	for (int i = 0; i < count - 1; i++)
	{
		DrawLine2D(vertices[i], vertices[i+1]);
	}
	DrawLine2D(vertices[count - 1], vertices[0]);
}

typedef struct _EdgeLUT
{
	Vector2 ver1;
	Vector2 ver2;
	float gradient;
	Colour4 colour;
} EdgeLUT;

void Rasterizer::ScanlineFillPolygon2D(const Vertex2d * vertices, int count)
{
	//TODO:
	//Ex 2.2 Implement the Rasterizer::ScanlineFillPolygon2D method method so that it is capable of drawing a solidly filled polygon.
	//Note: You can implement floodfill for this exercise however scanline fill is considered a more efficient and robust solution.
	//		You should be able to reuse DrawUnfilledPolygon2D here.
	//
	//Use Test 4 (Press F4) to test your solution, this is a simple test case as all polygons are convex.
	//Use Test 5 (Press F5) to test your solution, this is a complex test case with one non-convex polygon.

	//Ex 2.3 Extend Rasterizer::ScanlineFillPolygon2D method so that it is capable of alpha blending, i.e. draw translucent polygons.
	//Note: The variable mBlendMode indicates if the blend mode is set to alpha blending.
	//To do alpha blending during filling, the new colour of a point should be combined with the existing colour in the framebuffer using the alpha value.
	//Use Test 6 (Press F6) to test your solution
	
	int yMax = 0;
	int yMin = 1000;
	Colour4 colour = vertices[0].colour;
	std::vector<EdgeLUT> Edges;

	for (int i = 0; i < count - 1; i++)
	{
		EdgeLUT edge;
		if (vertices[i].position[1] > yMax)
			yMax = vertices[i].position[1];
		if (vertices[i].position[1] < yMin)
			yMin = vertices[i].position[1];

		if (yMax > mHeight)
			yMax = mHeight - 1;
		if (yMin < 1)
			yMin = 1;

		edge.gradient = ((float)vertices[i].position[0] - (float)vertices[i + 1].position[0]) / ((float)vertices[i].position[1] - (float)vertices[i + 1].position[1]);
		if (edge.gradient == -INFINITY)
			edge.gradient = 0;
		edge.ver1 = vertices[i].position;
		edge.ver2 = vertices[i + 1].position;
		Edges.push_back(edge);
	}
	EdgeLUT edge;
	edge.ver1 = vertices[count - 1].position;
	edge.ver2 = vertices[0].position;
	edge.gradient = ((float)edge.ver1[0] - (float)edge.ver2[0]) / ((float)edge.ver1[1] - (float)edge.ver2[1]);
	if (edge.gradient == -INFINITY)
	{
		edge.gradient = 0;
	}
	Edges.push_back(edge);

	for (int l = yMin; l < yMax; l++)
	{
		for (int e = 0; e < count; e++)
		{
			bool swapEdges = Edges[e].ver1[1] > Edges[e].ver2[1] ? true : false;
			if (swapEdges)
			{
				Vector2 temp(Edges[e].ver1);
				Edges[e].ver1 = Edges[e].ver2;
				Edges[e].ver2 = temp;
				
				if (Edges[e].gradient == -INFINITY)
				{
					Edges[e].gradient = 0;
				}
			}
			if ((Edges[e].ver1[1] <= l && l <= Edges[e].ver2[1]))
			{
				int x = Edges[e].ver1[0] + ((l - Edges[e].ver1[1]) * Edges[e].gradient); 

				bool duplicate = false;
				for (int i = 0; i < mScanlineLUT[l].size(); i++)
				{
					if (mScanlineLUT[l][i].pos_x == x)
					{
						if (Edges[e - 1].ver1[1] < Edges[e].ver1[1])
						{
							if (Edges[e].ver1[1] < Edges[e].ver2[1])
							{
								duplicate = true;
								break;
							}
						}
						else
						{
							if (Edges[e].ver1[1] > Edges[e].ver2[1])
							{
								duplicate = true;
								break;
							}
						}

					}

				}
				ScanlineLUTItem scanl = { vertices[0].colour, x };
				
				!duplicate ? mScanlineLUT[l].push_back(scanl) : duplicate = false;
			}
		}
		std::sort(mScanlineLUT[l].begin(), mScanlineLUT[l].end(), [](ScanlineLUTItem a, ScanlineLUTItem b) {return a.pos_x < b.pos_x; });
	}

	for (int scanline = yMin; scanline < yMax; scanline++)
	{
		int n = mScanlineLUT[scanline].size();
		if (n > 1)
		{
			for (int x = 0; x < n - 1; x+=2)
			{
				Vertex2d temp;
				temp.colour = mScanlineLUT[scanline][x].colour;
				temp.position = Vector2(mScanlineLUT[scanline][x].pos_x, scanline);
				Vertex2d temp2;
				temp2.colour = mScanlineLUT[scanline][x + 1].colour;
				temp2.position = Vector2(mScanlineLUT[scanline][x + 1].pos_x, scanline);

				DrawLine2D(temp, temp2);
			}
		}
	}
	ClearScanlineLUT();
}

void Rasterizer::ScanlineInterpolatedFillPolygon2D(const Vertex2d * vertices, int count)
{
	//TODO:
	//Ex 2.4 Implement Rasterizer::ScanlineInterpolatedFillPolygon2D method so that it is capable of performing interpolated filling.
	//Note: mFillMode is set to INTERPOLATED_FILL
	//		This exercise will be more straightfoward if Ex 1.3 has been implemented in DrawLine2D
	//Use Test 7 to test your solution
	int yMax = 0;
	int yMin = 1000;
	Colour4 colour = vertices[0].colour;
	std::vector<EdgeLUT> Edges;

	for (int i = 0; i < count - 1; i++)
	{
		EdgeLUT edge;
		if (vertices[i].position[1] > yMax)
			yMax = vertices[i].position[1];
		if (vertices[i].position[1] < yMin)
			yMin = vertices[i].position[1];

		if (yMax > mHeight)
			yMax = mHeight - 1;
		if (yMin < 1)
			yMin = 1;

		edge.gradient = ((float)vertices[i].position[0] - (float)vertices[i + 1].position[0]) / ((float)vertices[i].position[1] - (float)vertices[i + 1].position[1]);
		if (edge.gradient == -INFINITY)
			edge.gradient = 0;
		edge.ver1 = vertices[i].position;
		edge.ver2 = vertices[i + 1].position;
		edge.colour = vertices[i].colour;
		Edges.push_back(edge);
	}
	EdgeLUT edge;
	edge.ver1 = vertices[count - 1].position;
	edge.ver2 = vertices[0].position;
	edge.gradient = ((float)edge.ver1[0] - (float)edge.ver2[0]) / ((float)edge.ver1[1] - (float)edge.ver2[1]);
	edge.colour = vertices[count - 1].colour;
	if (edge.gradient == -INFINITY)
		edge.gradient = 0;
	Edges.push_back(edge);

	for (int l = yMin; l < yMax; l++)
	{
		for (int e = 0; e < count; e++)
		{
			bool swapEdges = Edges[e].ver1[1] > Edges[e].ver2[1] ? true : false;

			if (swapEdges)
			{
				Vector2 temp(Edges[e].ver1);
				Edges[e].ver1 = Edges[e].ver2;
				Edges[e].ver2 = temp;

				if (Edges[e].gradient == -INFINITY)
					Edges[e].gradient = 0;
			}

			if ((Edges[e].ver1[1] <= l && l <= Edges[e].ver2[1]))
			{
				int x = Edges[e].ver1[0] + ((l - Edges[e].ver1[1]) * Edges[e].gradient);

				bool duplicate = false;
				for (int i = 0; i < mScanlineLUT[l].size(); i++)
				{
					if (mScanlineLUT[l][i].pos_x == x)
					{
						if (Edges[e - 1].ver1[1] < Edges[e].ver1[1])
						{
							if (Edges[e].ver1[1] < Edges[e].ver2[1])
							{
								duplicate = true;
								break;
							}
						}
						else
						{
							if (Edges[e].ver1[1] > Edges[e].ver2[1])
							{
								duplicate = true;
								break;
							}
						}
					}
				}

				if (mFillMode == INTERPOLATED_FILLED)
				{
					Vector2 temp;
					temp[0] = l;
					temp[1] = x;

					Vertex2d temp2;
					temp2.position = Edges[e].ver1;
					Vertex2d temp3;
					temp3.position = Edges[e].ver2;

					float t = abs(temp2.position.Norm() - temp.Norm()) / abs(temp3.position.Norm() - temp2.position.Norm());

					if (e + 1 == count)
						colour = Edges[0].colour * t + Edges[e].colour * (1 - t);
					else
						colour = Edges[e+1].colour * t + Edges[e].colour * (1 - t);
				}
				ScanlineLUTItem scanl = { colour, x };

				!duplicate ? mScanlineLUT[l].push_back(scanl) : duplicate = false;
			}
		}

		std::sort(mScanlineLUT[l].begin(), mScanlineLUT[l].end(), [](ScanlineLUTItem a, ScanlineLUTItem b) {return a.pos_x < b.pos_x; });
	}

	for (int scanline = yMin; scanline < yMax; scanline++)
	{
		int n = mScanlineLUT[scanline].size();
		if (n > 1)
		{
			for (int x = 0; x < n - 1; x += 2)
			{
				Vertex2d temp;
				temp.colour = mScanlineLUT[scanline][x].colour;
				temp.position = Vector2(mScanlineLUT[scanline][x].pos_x, scanline);
				Vertex2d temp2;
				temp2.colour = mScanlineLUT[scanline][x + 1].colour;
				temp2.position = Vector2(mScanlineLUT[scanline][x + 1].pos_x, scanline);

				DrawLine2D(temp, temp2);
			}
		}
	}
	ClearScanlineLUT();
}

void Rasterizer::DrawCircle2D(const Circle2D & inCircle, bool filled)
{
	//TODO:
	//Ex 2.5 Implement Rasterizer::DrawCircle2D method so that it can draw a filled circle.
	//Note: For a simple solution, you can first attempt to draw an unfilled circle in the same way as drawing an unfilled polygon.
	//Use Test 8 to test your solution

	int nSegment = 24;
	const float PI = 3.1415926535897;

	int x = inCircle.centre[0];
	int y = inCircle.centre[1];
	int radius = inCircle.radius;
	float dt = 2 * PI / nSegment;

	Vertex2d temp;
	temp.position[0] = x + radius;
	temp.position[1] = y;
	temp.colour = inCircle.colour;
	Vertex2d temp2;
	temp2.colour = inCircle.colour;

	for (float t = 0; t < 2 * PI; t += dt)
	{
		temp2.position[0] = x + radius * cosf(t);
		temp2.position[1] = y + radius * sinf(t);
		DrawLine2D(temp, temp2);
		temp = temp2;
		dt = 2 * PI / nSegment;
	}
	temp2.position[0] = x + radius;
	temp2.position[1] = y;
	DrawLine2D(temp, temp2);
	/*temp2.position[0] = inCircle.centre[0];
	temp2.position[1] = inCircle.centre[1];
	DrawLine2D(temp, temp2);*/
}

Framebuffer *Rasterizer::GetFrameBuffer() const
{
	return mFramebuffer;
}
