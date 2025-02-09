/*
* Copyright(c) 2019, Youngsun Kwon, Donghyuk Kim, Inkyu An, and Sung-eui Yoon, KAIST
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met :
*
*     * Redistributions of source code must retain the above copyright notice, this
*       list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright notice,
*       this list of conditions and the following disclaimer in the documentation
*       and / or other materials provided with the distribution.
*     * Neither the name of SuperRay nor the names of its
*       contributors may be used to endorse or promote products derived from
*       this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

#ifndef GRIDMAP3D_CULLINGREGION_GRID3D_H
#define GRIDMAP3D_CULLINGREGION_GRID3D_H

#include <gridmap3D/gridmap3D.h>
#include <gridmap3D_superray/SuperRayGenerator.h>

#include <tr1/unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
# include <string.h>


using namespace std;
// using namespace chrono;

#define EPSS 1e-6
namespace gridmap3D{
    class CullingRegionGrid3D : public OccupancyGrid3DBase<Grid3DNode> {
    public:
        /// Default constructor, sets resolution of grid
        CullingRegionGrid3D(double resolution,double min_x,double max_x,double min_y,double max_y,double min_z,double max_z);

        /**
         * Reads a Grid3D from a binary file
         * @param _filename
         *
         */
        // CullingRegionGrid3D(std::string _filename);virtual 

        ~CullingRegionGrid3D(){

            // kown_boxfile.open("/home/jackykong/motionplanning/FUEL_ws/src/Exploration_sim/octomap_mapping/octomap_server/data/known_points.txt", std::ios_base::out);
            // int count = 0;
            // for (int i = 0; i<known_points.size(); ++i)
            // {
            //     kown_boxfile << known_points[i](0) << "	" << known_points[i](1) << " " << known_points[i](2) << endl;
            //     count++;
            // }
            // kown_boxfile << "end of file" << endl;
            // std::cout << "kown point write " << count <<std::endl;
            // kown_boxfile.close();
            myfile.close();
        };

        /// virtual constructor: creates a new object of same type
        /// (Covariant return type requires an up-to-date compiler)
        CullingRegionGrid3D* create() const { return new CullingRegionGrid3D(resolution,bbx_min,bbx_max,bby_min,bby_max,bbz_min,bbz_max); }

        std::string getGridType() const { return "CullingRegionGrid3D"; }

        // Super Rays and Culling Region based Updates

        /**
		 * Integrate a Pointcloud (in global reference frame), parallelized with OpenMP.
		 * This function simply inserts all rays of the point clouds, similar to insertPointCloudRays of gridmap3D::Grid3D.
		 * Occupied nodes have a preference over free ones.
		 *
		 * @param scan Pointcloud (measurement endpoints), in global reference frame
		 * @param origin measurement origin in global reference frame
		 */
        virtual void insertPointCloudRays(const Pointcloud& scan, const point3d& origin);

        /**
		 * Integrate a Pointcloud (in global reference frame) using SuperRay, parallelized with OpenMP.
		 * This function converts a point clouds into superrays, and then inserts all superrays out of the point clouds.
		 * Occupied nodes have a preference over free ones.
		 *
		 * @param scan Pointcloud (measurement endpoints), in global reference frame
		 * @param origin measurement origin in global reference frame
		 * @param threshold threshold for limiting to generate super rays
		 */
        virtual void insertSuperRayCloudRays(const Pointcloud& scan, const point3d& origin, const int threshold);

        virtual void trans_knownpoints(std::vector<point3d>& known_points_in)
        {
            known_points_in = known_points;
        }

    protected:
        /**
		 * Build a culling region by utilizing the occupancy information updated to the map.
		 * The implementation is based on a priority queue according to the Manhattan distance from the origin cell.
         *
		 * @param origin measurement origin in global reference frame
		 * @param max_propagation maximum level of propagation; Manhattan distance from the origin cell
         * @return culling region limited by the maximum level of the propagation
		 */
        KeySet buildCullingRegion(const point3d& origin, const int max_propagation);

        /**
		 * Build a culling region by utilizing the occupancy information updated to the map.
         *
         * @param scan Pointcloud (measurement endpoints), in global reference frame
		 * @param origin measurement origin in global reference frame
         * @return culling region limited by the range of measurements
		 */
        KeySet buildCullingRegion(const Pointcloud& scan, const point3d& origin);

        /**
		 * Build a culling region by utilizing the occupancy information updated to the map.
         *
         * @param superrays Super rays computed from the measurements in global reference frame
		 * @param origin measurement origin in global reference frame
         * @return culling region limited by the range of measurements
		 */
        KeySet buildCullingRegion(const SuperRayCloud& superrays, const point3d& origin);

        /**
         * Traces a sensor ray from origin (excluding) to end in the inverse direction (see the description),
         * returning Grid3DKeys of all nodes traversed by the beam.
         * During the traversal, the culling region stops the traversal when the ray
         * encounters the region for reducing the unnecessary traversals and updates.
         *
         * @param origin start coordinate of ray (end point of sensor ray)
         * @param end end coordinate of ray (sensor origin)
         * @param ray KeyRay structure that holds the keys of all nodes traversed by the ray, excluding the origin cell
         * @return Success of operation. Returning false usually means that one of the coordinates is out of the Grid3D's range
         */
        bool computeInverseRayKeys(const point3d& origin, const point3d& end, KeyRay& ray, KeySet& cullingregion);


        /**
         * Static member object which ensures that this Grid3D's prototype
         * ends up in the classIDMapping only once. You need this as a
         * static member in any derived grid3D class in order to read .og3
         * files through the AbstractGrid3D factory. You should also call
         * ensureLinking() once from the constructor.
         */
        class StaticMemberInitializer{
        public:
            StaticMemberInitializer() {
                CullingRegionGrid3D* grid = new CullingRegionGrid3D(0.1,0,10,0,10,0,10);
                grid->clearKeyRays();
                AbstractGrid3D::registerGridType(grid);
            }

            /**
             * Dummy function to ensure that MSVC does not drop the
             * StaticMemberInitializer, causing this grid failing to register.
             * Needs to be called from the constructor of this grid3D.
             */
            void ensureLinking() {};
        };
        /// static member to ensure static initialization (only once)
        static StaticMemberInitializer cullingregionGrid3DMemberInit;

    private:
        double m_res;
        int explored_voxelcount = 0;
        int curexpl_voxelcount = 0;
        std::tr1::unordered_map<int, int> point_hashmap;
        std::vector<point3d> known_points;
        int known_points_count = 0;
        std::string pkg_path;
        std::ofstream myfile, kown_boxfile;
        std::chrono::_V2::system_clock::time_point start;

        double bbx_min,bbx_max,bby_min,bby_max,bbz_min,bbz_max;

        // bool*** gridmap_mat = NULL;
        // int cube_numx,cube_numy,cube_numz;
    };
}

#endif
