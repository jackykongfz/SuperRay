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

#include <gridmap3D_cullingregion/CullingRegionGrid3D.h>

namespace gridmap3D{
    CullingRegionGrid3D::CullingRegionGrid3D(double in_resolution, double min_x,double max_x,double min_y,double max_y,double min_z,double max_z)
            : OccupancyGrid3DBase<Grid3DNode>(in_resolution) {
        cullingregionGrid3DMemberInit.ensureLinking();
        m_res = in_resolution;
        myfile.open("/home/jackykong/motionplanning/FUEL_ws/src/Exploration_sim/octomap_mapping/octomap_server/data/superray_data.txt", std::ios_base::out);//, std::ios_base::out
        start = std::chrono::system_clock::now();
        bbx_min = min_x;
        bbx_max = max_x;
        bby_min = min_y;
        bby_max = max_y;
        bbz_min = min_z;
        bbz_max = max_z;
    };

    /*CullingRegionGrid3D::CullingRegionGrid3D(std::string _filename)
            : OccupancyGrid3DBase<Grid3DNode>(0.1)  { // resolution will be set according to grid file
        readBinary(_filename);
    }*/

    CullingRegionGrid3D::StaticMemberInitializer CullingRegionGrid3D::cullingregionGrid3DMemberInit;

    void CullingRegionGrid3D::insertPointCloudRays(const Pointcloud& pc, const point3d& origin)
    {
        if (pc.size() < 1)
            return;

        // Build a culling region
        KeySet cullingregion = buildCullingRegion(pc, origin);

        // std::cout<<"origin = "<<origin(0)<<", " << origin(1) <<", "<<origin(2)<<std::endl;
        // std::cout<<"bounding box area = "<<bbx_min<<", " << bbx_max <<", "<<bby_min 
        //            << ", " << bby_max << ", " << bbz_min << ", " << bbz_max<<std::endl;

            double hash_cubesize = m_res;
            int cube_numx,cube_numy,cube_numz;
            cube_numx = 500/hash_cubesize;
            cube_numy = 500/hash_cubesize;
            cube_numz = 100/hash_cubesize;
            curexpl_voxelcount = 0;

            

        // Update the occupancies of the map
#ifdef _OPENMP
        omp_set_num_threads(this->keyrays.size());
#pragma omp parallel for
#endif
        for (int i = 0; i < (int)pc.size(); ++i) {
            const point3d& p = pc[i];
            unsigned threadIdx = 0;
#ifdef _OPENMP
            threadIdx = omp_get_thread_num();
#endif
            KeyRay* keyray = &(this->keyrays.at(threadIdx));


            if (this->computeInverseRayKeys(p, origin, *keyray, cullingregion)){
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    // Update the traversed cells to have the free states
                    for (KeyRay::iterator it = keyray->begin(); it != keyray->end(); ++it){
                        updateNode(*it, false);
                        point3d keys_pt = keyToCoord(*it);

                        // std::cout<<"current pt = "<<keys_pt(0)<<", " << keys_pt(1) <<", "<<keys_pt(2)<<std::endl;

                        //add bounding box 
                        if(keys_pt(0) > bbx_max || keys_pt(0) < bbx_min)
                        {
                            continue;
                        }else if(keys_pt(1) > bby_max || keys_pt(1) < bby_min)
                        {
                            continue;
                        }else if(keys_pt(2) > bbz_max || keys_pt(2) < bbz_min)
                        {
                            continue;
                        }

                        //check if hashmap has value
                        int ind_x = (round((keys_pt(0) - (bbx_min - 50) + EPSS)/hash_cubesize));
                        int ind_y = (round((keys_pt(1) - (bby_min - 50) + EPSS)/hash_cubesize));
                        int ind_z = (round((keys_pt(2) - (bbz_min - 50) + EPSS)/hash_cubesize));

                        if (ind_x < 0 || ind_y < 0 || ind_z <0)
                        {
                            std::cout << "hash value = " << ind_x << ", " << ind_y << ", " << ind_z <<std::endl;
                        }

                        long int box_index = ind_x + ind_y*cube_numx + ind_z*cube_numx*cube_numy;
                        // point_hashmap.insert(pair<int, >(box_index, pt_in));
                        if(point_hashmap[box_index] > 0)
                        {
                            continue;
                        }else{
                            point_hashmap[box_index] = 1;
                            curexpl_voxelcount++;

                            point3d center_pt;
                            center_pt(0) = ind_x * hash_cubesize + (bbx_min - 50);
                            center_pt(1) = ind_y * hash_cubesize + (bby_min - 50);
                            center_pt(2) = ind_z * hash_cubesize + (bbz_min - 50);

                            known_points.push_back(center_pt);
                        }
                    }
                }
            }
        }

        // Update the cells containing the end points to have the occupied states
        for (int i = 0; i < (int)pc.size(); ++i){
            updateNode(pc[i], true);
                        point3d keys_pt = pc[i];

                        //add bounding box 
                        if(keys_pt(0) > bbx_max || keys_pt(0) < bbx_min)
                        {
                            continue;
                        }else if(keys_pt(1) > bby_max || keys_pt(1) < bby_min)
                        {
                            continue;
                        }else if(keys_pt(2) > bbz_max || keys_pt(2) < bbz_min)
                        {
                            continue;
                        }

                        //check if hashmap has value
                        int ind_x = (round((keys_pt(0) - (bbx_min - 50) + EPSS)/hash_cubesize));
                        int ind_y = (round((keys_pt(1) - (bby_min - 50) + EPSS)/hash_cubesize));
                        int ind_z = (round((keys_pt(2) - (bbz_min - 50) + EPSS)/hash_cubesize));

                        if (ind_x < 0 || ind_y < 0 || ind_z <0)
                        {
                            std::cout << "hash value = " << ind_x << ", " << ind_y << ", " << ind_z <<std::endl;
                        }

                        long int box_index = ind_x + ind_y*cube_numx + ind_z*cube_numx*cube_numy;
                        // point_hashmap.insert(pair<int, >(box_index, pt_in));
                        if(point_hashmap[box_index] > 0)
                        {
                            continue;
                        }else{
                            point_hashmap[box_index] = 1;
                            curexpl_voxelcount++;

                            point3d center_pt;
                            center_pt(0) = ind_x * hash_cubesize + (bbx_min - 50);
                            center_pt(1) = ind_y * hash_cubesize + (bby_min - 50);
                            center_pt(2) = ind_z * hash_cubesize + (bbz_min - 50);

                            known_points.push_back(center_pt);
                        }
        }

        explored_voxelcount = explored_voxelcount + curexpl_voxelcount;
        std::cout << "Explored voxel count = " <<  explored_voxelcount << std::endl;
        auto end2 = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start);
        myfile << double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den  << " " << curexpl_voxelcount << " " << explored_voxelcount << " " << m_res*m_res*m_res*explored_voxelcount << std::endl;

    }

    void CullingRegionGrid3D::insertSuperRayCloudRays(const Pointcloud& pc, const point3d& origin, const int threshold)
    {
        if (pc.size() < 1)
            return;

        // Generate the super rays
        SuperRayGenerator srgenerator(resolution, grid_max_val, threshold);
        SuperRayCloud srcloud;
        srgenerator.GenerateSuperRay(pc, origin, srcloud);

        // Build a culling region
        KeySet cullingregion = buildCullingRegion(srcloud, origin);

        // Update the occupancies of the map
#ifdef _OPENMP
        omp_set_num_threads(this->keyrays.size());
#pragma omp parallel for
#endif
        for (int i = 0; i < (int)srcloud.size(); ++i) {
            const point3d& p = srcloud[i].p;
            unsigned threadIdx = 0;
#ifdef _OPENMP
            threadIdx = omp_get_thread_num();
#endif
            KeyRay* keyray = &(this->keyrays.at(threadIdx));

            if (this->computeInverseRayKeys(p, origin, *keyray, cullingregion)){
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    // Update the traversed cells to have the free states
                    for (KeyRay::iterator it = keyray->begin(); it != keyray->end(); it++) {
                        updateNode(*it, prob_miss_log * srcloud[i].w);
                    }
                }
            }
        }

        // Update the cells containing the end points to have the occupied states
        for (int i = 0; i < (int)srcloud.size(); ++i){
            updateNode(srcloud[i].p, prob_hit_log * srcloud[i].w);
        }
    }

    KeySet CullingRegionGrid3D::buildCullingRegion(const point3d& origin, const int max_propagation)
    {
        KeySet cullingregion;

        Grid3DKey originKey = coordToKey(origin);

        KeySet* cur_candidates = new KeySet;
        cur_candidates->insert(originKey);

        for(int cur_level = 0; cur_level <= max_propagation; cur_level++){
            KeySet* next_candidates = new KeySet;
            for(KeySet::iterator it = cur_candidates->begin(); it != cur_candidates->end(); it++){
                // Check the insertion of the cell into the culling region
                const Grid3DKey& key = *it;

                // The first condition: does the cell have a fully free state?
                int step[3] = {0, 0, 0};
                Grid3DNode* node = search(key);
                if (!node || node->getLogOdds() > clamping_thres_min)
                    continue;

                // The second condition: are all the neighbor cells in the culling region?
                bool insertion = true;
                for (int axis = 0; axis < 3; axis++){
                    // Find a neighbor cell in the direction of the axis
                    if (key[axis] > originKey[axis])		step[axis] = -1;
                    else if (key[axis] < originKey[axis])	step[axis] = 1;

                    if (step[axis] != 0){
                        // Check a neighbor cell
                        Grid3DKey checkKey = key;
                        checkKey[axis] += step[axis];

                        // The neighbor cell is not in culling region
                        if (cullingregion.find(checkKey) == cullingregion.end()){
                            insertion = false;
                            break;
                        }
                    }
                }

                // Insert the cell into the culling region
                if(!insertion)
                    continue;
                cullingregion.insert(key);

                // Find the candidates in the next level
                if(cur_level != max_propagation){
                    for (int axis = 0; axis < 3; axis++){
                        Grid3DKey candidate = key;

                        if (step[axis] != 0){
                            candidate[axis] += (-step[axis]);
                            next_candidates->insert(candidate);
                        }
                        else{
                            candidate[axis] += 1;
                            next_candidates->insert(candidate);

                            candidate[axis] += -2;
                            next_candidates->insert(candidate);
                        }
                    }
                }
            }

            // Propagation: move to the next level
            delete cur_candidates;
            cur_candidates = next_candidates;
            if(cur_candidates->size() <= 0)
                break;
        }

        delete cur_candidates;

        return cullingregion;
    }

    KeySet CullingRegionGrid3D::buildCullingRegion(const Pointcloud& pc, const point3d& origin)
    {
        // Find the maximum distance between the sensor origin and the end point
        double max_dist = 0.0;
        for(int i = 0; i < (int)pc.size(); i++){
            const point3d& p = pc[i];
            double distance = (p - origin).norm();
            if(distance > max_dist)
                max_dist = distance;
        }

        // Build a culling region limited the minimum distance
        return buildCullingRegion(origin, (int)(max_dist / resolution));
    }

    KeySet CullingRegionGrid3D::buildCullingRegion(const SuperRayCloud& superrays, const point3d& origin)
    {
        // Find the maximum distance between the sensor origin and the end point
        double max_dist = 0.0;
        for(int i = 0; i < (int)superrays.size(); i++){
            const point3d& p = superrays[i].p;
            double distance = (p - origin).norm();
            if(distance > max_dist)
                max_dist = distance;
        }

        // Build a culling region limited the minimum distance
        return buildCullingRegion(origin, (int)(max_dist / resolution));
    }

    bool CullingRegionGrid3D::computeInverseRayKeys(const point3d& origin, const point3d& end, KeyRay& ray, KeySet& cullingregion)
    {
        // see "A Faster Voxel Traversal Algorithm for Ray Tracing" by Amanatides & Woo
        // basically: DDA in 3D

        ray.reset();

        Grid3DKey key_origin, key_end;
        if ( !coordToKeyChecked(origin, key_origin) || !coordToKeyChecked(end, key_end) ) {
            GRIDMAP3D_WARNING_STR("coordinates ( " << origin << " -> " << end << ") out of bounds in computeRayKeys");
            return false;
        }

        if (key_origin == key_end)
            return true; // same tree cell, we're done.

        // Initialization phase -------------------------------------------------------
        point3d direction = (end - origin);
        float length = (float) direction.norm();
        direction /= length; // normalize vector

        int    step[3];
        double tMax[3];
        double tDelta[3];

        Grid3DKey current_key = key_origin;

        for(unsigned int i = 0; i < 3; ++i) {
            // compute step direction
            if (direction(i) > 0.0)         step[i] = 1;
            else if (direction(i) < 0.0)    step[i] = -1;
            else                            step[i] = 0;

            // compute tMax, tDelta
            if (step[i] != 0) {
                // corner point of voxel (in direction of ray)
                double voxelBorder = this->keyToCoord(current_key[i]);
                voxelBorder += (float) (step[i] * this->resolution * 0.5);

                tMax[i] = ( voxelBorder - origin(i) ) / direction(i);
                tDelta[i] = this->resolution / fabs( direction(i) );
            }
            else {
                tMax[i] =  std::numeric_limits<double>::max( );
                tDelta[i] = std::numeric_limits<double>::max( );
            }
        }

        // Incremental phase  ---------------------------------------------------------
        bool done = false;
        while (!done) {
            // find minimum tMax:
            unsigned int dim;
            if (tMax[0] < tMax[1]){
                if (tMax[0] < tMax[2]) dim = 0;
                else                   dim = 2;
            }
            else {
                if (tMax[1] < tMax[2]) dim = 1;
                else                   dim = 2;
            }

            // advance in direction "dim"
            current_key[dim] += step[dim];
            tMax[dim] += tDelta[dim];

            assert (current_key[dim] < 2*this->grid_max_val);

            // Culling out the traversal
            if (cullingregion.find(current_key) != cullingregion.end()){
                done = true;
                break;
            }
            // reached endpoint, key equv?
            else if (current_key == key_end) {
                ray.addKey(current_key);
                done = true;
                break;
            }
            else {
                // reached endpoint world coords?
                // dist_from_origin now contains the length of the ray when traveled until the border of the current voxel
                double dist_from_origin = std::min(std::min(tMax[0], tMax[1]), tMax[2]);
                // if this is longer than the expected ray length, we should have already hit the voxel containing the end point with the code above (key_end).
                // However, we did not hit it due to accumulating discretization errors, so this is the point here to stop the ray as we would never reach the voxel key_end
                if (dist_from_origin > length) {
                    done = true;
                    break;
                }
                else {  // continue to add freespace cells
                    ray.addKey(current_key);
                }
            }

            assert ( ray.size() < ray.sizeMax() - 1);
        } // end while

        return true;
    }
}