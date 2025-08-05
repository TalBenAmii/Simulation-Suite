#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <math.h>
#include <stdint.h>
#include <cmath> // Include the cmath library for pow
using namespace std;

using std::FILE;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ifstream;
using std::stringstream;


class Block;
class Cache;


typedef struct {
	unsigned instruction_count;
	unsigned L1_hit_count;
	unsigned L2_hit_count;
    unsigned long long t_access;
} SIM_stats;

class Block
{
    friend class Cache;
public:
    unsigned addr;
    unsigned tag;
    bool dirty;
    bool invalid;
    unsigned LRU;
    
    Block(unsigned addr = 0, unsigned tag = 0, unsigned LRU = 0)
        : addr(addr), tag(tag), dirty(false), invalid(true), LRU(LRU) 
    {}
};


// cache constants
const static int BLOCK_BITS_SIZE = 32;
const static int ALIGNED_BITS = 2;
const static int MISS = -1;
const static int WRITE_ALOOCATE = 1;
const static int WRITE_NO_ALOOCATE = 0;
class Cache
{
    friend class Block;
public:
    
    unsigned    BLOCK_SIZE,
                OFFSET_BITS,
                L1_NUM_SETS,
                L1_NUM_WAYS,
                L1_CACHE_SIZE,
                L1_NUM_BLOCKS,
                L1_SET_BITS,
                L1_TAG_BITS,  
                L2_NUM_SETS,
                L2_NUM_WAYS,
                L2_CACHE_SIZE,
                L2_NUM_BLOCKS,
                L2_SET_BITS,
                L2_TAG_BITS;
    
    // attributes
    unsigned    MemCyc,
                BSize,
                L1Size,
                L2Size,
                L1Assoc,
                L2Assoc, 
                L1Cyc,
                L2Cyc,
                WrAlloc;

    SIM_stats stats;
    vector< vector<Block> > L1, L2;

    

    Cache() = default;
    ~Cache() = default;
    Cache(  unsigned MemCyc, unsigned BSize, unsigned L1Size, unsigned L2Size,
            unsigned L1Assoc, unsigned L2Assoc, unsigned L1Cyc, unsigned L2Cyc,
            unsigned WrAlloc);

    void init_LRUs();
    int find_way_in_level(unsigned addr, int level);
    bool level_full(int level, int set);
    void dirty_block(int level, int way, int set);
    void invalidate_block(int level, int way, int set);
    void update_LRU(int level, int way, int set);
    void update_lower_level(unsigned addr);
    int choose_victim_way(int level, int set);
    
    void add_to_level(unsigned addr, int level);
    void add_to_cache(unsigned addr);
    void Read_L1_hit(unsigned addr, int way_L1);
    void Read_L1_miss_L2_hit(unsigned addr, int way_L2);
    void Read_L1_miss_L2_miss(unsigned addr);
    void Write_L1_hit(unsigned addr, int way_L1);
    void Write_L1_miss_L2_hit(unsigned addr, int way_L2);
    void Write_L1_miss_L2_miss(unsigned addr);

    void print_cache(int iteration);
    void Read(unsigned addr);
    void Write(unsigned addr);
    double get_L1_miss_rate();
    double get_L2_miss_rate();
    double get_avg_access_time();
};

/**
 * @brief Constructs a Cache object with the specified parameters.
 *
 * @param MemCyc Memory cycle time.
 * @param BSize Block size in bits.
 * @param L1Size L1 cache size in bits.
 * @param L2Size L2 cache size in bits.
 * @param L1Assoc L1 cache associativity in bits.
 * @param L2Assoc L2 cache associativity in bits.
 * @param L1Cyc L1 cache cycle time.
 * @param L2Cyc L2 cache cycle time.
 * @param WrAlloc Write allocate policy (1 for write-allocate, 0 for no write-allocate).
 */
Cache::Cache(unsigned MemCyc, unsigned BSize, unsigned L1Size, unsigned L2Size,
             unsigned L1Assoc, unsigned L2Assoc, unsigned L1Cyc, unsigned L2Cyc,
             unsigned WrAlloc)
    : MemCyc(MemCyc), BSize(BSize), L1Size(L1Size),
      L2Size(L2Size), L1Assoc(L1Assoc), L2Assoc(L2Assoc),
      L1Cyc(L1Cyc), L2Cyc(L2Cyc), WrAlloc(WrAlloc)
{
    stats = {0, 0, 0, 0};
    
    
    BLOCK_SIZE = pow(2, BSize);
    OFFSET_BITS = BSize;
    L1_NUM_SETS = pow(2, L1Size - BSize - L1Assoc);
    L1_NUM_WAYS = pow(2, L1Assoc);
    L1_CACHE_SIZE = pow(2, L1Size);
    L1_NUM_BLOCKS = L1_NUM_SETS * L1_NUM_WAYS;
    L1_SET_BITS = L1Size - BSize - L1Assoc;
    L1_TAG_BITS = BLOCK_BITS_SIZE - OFFSET_BITS - L1_SET_BITS;    
    L2_NUM_SETS = pow(2, L2Size - BSize - L2Assoc);
    L2_NUM_WAYS = pow(2, L2Assoc);
    L2_CACHE_SIZE = pow(2, L2Size);
    L2_NUM_BLOCKS = L2_NUM_SETS * L2_NUM_WAYS;
    L2_SET_BITS = L2Size - BSize - L2Assoc;
    L2_TAG_BITS = BLOCK_BITS_SIZE - OFFSET_BITS - L2_SET_BITS;

    L1.resize(L1_NUM_WAYS, vector<Block>(L1_NUM_SETS));
    L2.resize(L2_NUM_WAYS, vector<Block>(L2_NUM_SETS));
    init_LRUs();
}

/**
 * @brief Initializes the Least Recently Used (LRU) values and invalidates all cache lines
 *        in both L1 and L2 caches.
 * 
 * This function sets up the initial state of the cache by assigning an LRU value to each
 * cache line and marking all cache lines as invalid. The LRU values are assigned in a 
 * sequential manner starting from 0 for each way in both L1 and L2 caches.
 * 
 * The function iterates over all ways and sets in the L1 cache, assigning an LRU value 
 * and setting the invalid flag to true for each cache line. It then repeats the process 
 * for the L2 cache.
 * 
 * @note The constants L1_NUM_WAYS, L1_NUM_SETS, L2_NUM_WAYS, and L2_NUM_SETS are assumed 
 *       to be defined elsewhere in the code.
 */
void Cache::init_LRUs()
{
    unsigned curr_LRU = 0;
    for (int i = 0; i < L1_NUM_WAYS; i++) {
        for (int j = 0; j < L1_NUM_SETS; j++) {
            L1[i][j].LRU = curr_LRU;
            L1[i][j].invalid = true;
        }
        curr_LRU++;
    }
    curr_LRU = 0;
    for (int i = 0; i < L2_NUM_WAYS; i++) {
        for (int j = 0; j < L2_NUM_SETS; j++) {
            L2[i][j].LRU = curr_LRU;
            L2[i][j].invalid = true;
        }
        curr_LRU++;
    }
}

/**
 * @brief Finds the way in the specified cache level that contains the given address.
 *
 * This function searches through the specified cache level (L1 or L2) to find the way
 * that contains the block with the given address. It checks each way in the set
 * corresponding to the address and returns the index of the way if a match is found.
 *
 * @param addr The address to search for in the cache.
 * @param level The cache level to search in (1 for L1, 2 for L2).
 * @return The index of the way containing the address if found, otherwise returns MISS.
 */
int Cache::find_way_in_level(unsigned addr, int level)
{
    vector<vector<Block> > &curr_level = (level == 1) ? L1 : L2;
    int curr_level_ways = (level == 1) ? L1_NUM_WAYS : L2_NUM_WAYS;
    int curr_set = (addr >> (OFFSET_BITS)) % (level == 1 ? L1_NUM_SETS : L2_NUM_SETS);
    unsigned curr_tag = (addr >> (OFFSET_BITS + (level == 1 ? L1_SET_BITS : L2_SET_BITS)));
    for (int i = 0; i < curr_level_ways; i++)
    {
        if (curr_level[i][curr_set].tag == curr_tag && !(curr_level[i][curr_set].invalid))
        {
            return i;
        }
    }
    return MISS;
}

/**
 * @brief Checks if a given cache level and set is full.
 *
 * This function determines whether all the blocks in a specified set of a given cache level
 * are valid (i.e., not invalid). It iterates through all the ways in the specified set and
 * checks the 'invalid' flag of each block.
 *
 * @param level The cache level to check (1 for L1, 2 for L2).
 * @param set The set index within the specified cache level to check.
 * @return true if all blocks in the specified set are valid (i.e., the set is full), false otherwise.
 */
bool Cache::level_full(int level, int set)
{
    vector<vector<Block> > &curr_level = (level == 1) ? L1 : L2;
    int curr_level_ways = (level == 1) ? L1_NUM_WAYS : L2_NUM_WAYS;
    for (int i = 0; i < curr_level_ways; i++)
    {
        if (curr_level[i][set].invalid)
        {
            return false;
        }
    }
    return true;
}

/**
 * Marks a block in the cache as dirty.
 *
 * @param level The cache level (1 for L1, 2 for L2).
 * @param way The way index in the cache set.
 * @param set The set index in the cache.
 */
void Cache::dirty_block(int level, int way, int set)
{
    vector<vector<Block> > &curr_level = (level == 1) ? L1 : L2;
    curr_level[way][set].dirty = true;
}

/**
 * @brief Invalidates a cache block in the specified cache level.
 *
 * This function marks a cache block as invalid in either the L1 or L2 cache,
 * based on the provided cache level, way, and set indices.
 *
 * @param level The cache level to invalidate the block in (1 for L1, 2 for L2).
 * @param way The way index of the cache block to invalidate.
 * @param set The set index of the cache block to invalidate.
 */
void Cache::invalidate_block(int level, int way, int set)
{
    vector<vector<Block> > &curr_level = (level == 1) ? L1 : L2;
    curr_level[way][set].invalid = true;
}

/**
 * @brief Updates the Least Recently Used (LRU) status of a cache block.
 *
 * This function updates the LRU status of a cache block in either the L1 or L2 cache level.
 * It sets the LRU value of the specified block to the highest value (most recently used)
 * and decrements the LRU values of other blocks in the same set that have higher LRU values.
 *
 * @param level The cache level to update (1 for L1, 2 for L2).
 * @param way The way index of the cache block to update.
 * @param set The set index of the cache block to update.
 */
void Cache::update_LRU(int level, int way, int set)
{
    vector<vector<Block> > &curr_level = (level == 1) ? L1 : L2;
    int curr_level_ways = (level == 1) ? L1_NUM_WAYS : L2_NUM_WAYS;
    int x = curr_level[way][set].LRU;
    curr_level[way][set].LRU = curr_level_ways - 1;
    for (int i = 0; i < curr_level_ways; i++)
    {
        if (i != way && curr_level[i][set].LRU > x)
        {
            curr_level[i][set].LRU--;
        }
    }
}

/**
 * @brief Updates the lower level cache with the given address.
 *
 * This function updates the lower level cache (L2) with the specified address.
 * It calculates the set index and way index for the address in the L2 cache.
 * If the address is not found in the L2 cache (miss), the function returns.
 * Otherwise, it marks the block as dirty and updates the LRU (Least Recently Used) status.
 *
 * @param addr The address to be updated in the lower level cache.
 */
void Cache::update_lower_level(unsigned addr)
{
    int curr_set = (addr >> OFFSET_BITS) % L2_NUM_SETS;
    int curr_way = find_way_in_level(addr, 2);
    if (curr_way == MISS)
    {
        return;
    }
    dirty_block(2, curr_way, curr_set);
    update_LRU(2, curr_way, curr_set);    
}

/**
 * @brief Chooses a victim way in the cache set for replacement based on the LRU (Least Recently Used) policy.
 *
 * This function iterates through the cache ways in the specified level and set, and selects the way
 * that has the LRU value of 0, indicating it is the least recently used and should be replaced.
 *
 * @param level The cache level (1 for L1 cache, 2 for L2 cache).
 * @param set The index of the cache set to search for a victim way.
 * @return The index of the chosen victim way, or MISS if no victim way is found.
 */
int Cache::choose_victim_way(int level, int set)
{
    int chosen_way = MISS;
    vector<vector<Block> > &curr_level = (level == 1) ? L1 : L2;
    int curr_level_ways = (level == 1) ? L1_NUM_WAYS : L2_NUM_WAYS;
    for (int i = 0; i < curr_level_ways; i++)
    {
        if (curr_level[i][set].LRU == 0)
        {
            chosen_way = i;
            return chosen_way;
        }
    }
    return chosen_way;
}

/**
 * @brief Reads data from the cache given an address.
 *
 * This function increments the instruction count and attempts to read data from the cache
 * at the specified address. It first checks for a hit in the L1 cache. If a hit is found,
 * it processes the L1 hit. If there is a miss in the L1 cache but a hit in the L2 cache,
 * it processes the L1 miss and L2 hit. If there is a miss in both the L1 and L2 caches,
 * it processes the L1 and L2 miss.
 *
 * @param addr The address from which to read data.
 */
void Cache::Read(unsigned addr)
{
    stats.instruction_count++;
    int way_L1 = find_way_in_level(addr, 1);
    int way_L2 = find_way_in_level(addr, 2);

    if (way_L1 != MISS)  // L1 hit
    {
        Cache::Read_L1_hit(addr, way_L1);
        return;
    }
    else if (way_L2 != MISS)  // L2 hit && L1 miss
    {
        Cache::Read_L1_miss_L2_hit(addr, way_L2);
        return;
    }
    else  // L2 miss && L1 miss
    {
        Cache::Read_L1_miss_L2_miss(addr);
        return;
    }    
}

/**
 * @brief Handles a read hit in the L1 cache.
 *
 * This function updates the cache statistics and the LRU (Least Recently Used) 
 * information when a read hit occurs in the L1 cache.
 *
 * @param addr The address of the memory being accessed.
 * @param way_L1 The way index in the L1 cache where the hit occurred.
 */
void Cache::Read_L1_hit(unsigned addr, int way_L1)
{
    // update stats
    stats.L1_hit_count++;
    stats.t_access += L1Cyc;
    
    int set_L1 = (addr >> OFFSET_BITS) % L1_NUM_SETS;
    update_LRU(1, way_L1, set_L1);
    return;
}

/**
 * @brief Handles the scenario where a read operation results in an L1 cache miss but an L2 cache hit.
 * 
 * This function updates the cache statistics, manages the LRU (Least Recently Used) policy,
 * and handles the eviction and addition of cache blocks in the L1 cache.
 * 
 * @param addr The address of the memory being accessed.
 * @param way_L2 The way index in the L2 cache where the block was found.
 */
void Cache::Read_L1_miss_L2_hit(unsigned addr, int way_L2)
{
    // update stats
    stats.L2_hit_count++;
    stats.t_access += L1Cyc;
    stats.t_access += L2Cyc;
    
    int set_L1 = (addr >> OFFSET_BITS) % L1_NUM_SETS;
    int set_L2 = (addr >> OFFSET_BITS) % L2_NUM_SETS;
    update_LRU(2, way_L2, set_L2);
    
    if (level_full(1, set_L1))
    {
        int victim_way = choose_victim_way(1, set_L1);
        if (L1[victim_way][set_L1].dirty)
        {
            update_lower_level(L1[victim_way][set_L1].addr);
        }
        invalidate_block(1, victim_way, set_L1);
    }
    add_to_level(addr, 1); // supposed to add only to L1 and to update LRU
    return;
}

/**
 * @brief Handles a cache read operation when both L1 and L2 caches miss.
 *
 * This function updates the access time statistics to account for the time
 * taken to access L1 cache, L2 cache, and main memory. It then adds the
 * requested address to the cache.
 *
 * @param addr The memory address that caused the cache miss.
 */
void Cache::Read_L1_miss_L2_miss(unsigned addr)
{
    // update stats
    stats.t_access += L1Cyc;
    stats.t_access += L2Cyc;
    stats.t_access += MemCyc;

    add_to_cache(addr);
    return;
}

/**
 * @brief Adds a block to the cache.
 *
 * This function adds a block to both L1 and L2 caches. If the cache set is full,
 * it will choose a victim block to evict based on the replacement policy.
 * If the victim block is dirty, it will update the lower level cache or memory.
 * It also ensures coherence by invalidating the corresponding block in the other cache level.
 *
 * @param addr The address of the block to be added to the cache.
 */
void Cache::add_to_cache(unsigned addr)
{
    int set_L1 = (addr >> OFFSET_BITS) % L1_NUM_SETS;
    int set_L2 = (addr >> OFFSET_BITS) % L2_NUM_SETS;

    // add to L2
    if(level_full(2, set_L2))
    {
        int victim_way = choose_victim_way(2, set_L2);
        Block& victim_block = L2[victim_way][set_L2];
        // if (victim_block.dirty){ stats.t_access += MemCyc;} // write back to memory (not sure if needed)
        invalidate_block(2, victim_way, set_L2);

        // snoop above (L1)
        int invalid_block_set_L1 = (victim_block.addr >> OFFSET_BITS) % L1_NUM_SETS;
        int invalid_block_way_L1 = find_way_in_level(victim_block.addr, 1);
        if (invalid_block_way_L1 != MISS) // victim found in L1
        {
            // if (invalid_block.dirty){ stats.t_access += MemCyc;} // write back to memory (not sure if needed)
            invalidate_block(1, invalid_block_way_L1, invalid_block_set_L1);
        }
    }
    add_to_level(addr, 2);

    // add to L1
    if(level_full(1, set_L1))
    {
        int victim_way = choose_victim_way(1, set_L1);
        Block& victim_block = L1[victim_way][set_L1];

        if (victim_block.dirty)
        {
            update_lower_level(victim_block.addr);
        }
        invalidate_block(1, victim_way, set_L1);
    }
    add_to_level(addr, 1);
}

/**
 * @brief Adds a block to the specified cache level.
 *
 * This function attempts to add a block with the given address to the specified cache level (L1 or L2).
 * It first determines the appropriate set and tag for the address, then searches for an invalid block
 * within the set to store the new block. If an invalid block is found, it updates the block's address,
 * tag, and status, and updates the LRU (Least Recently Used) information.
 *
 * @param addr The address of the block to be added to the cache.
 * @param level The cache level to which the block should be added (1 for L1, 2 for L2).
 */
void Cache::add_to_level(unsigned addr, int level)
{
    vector<vector<Block> > &curr_level = (level == 1) ? L1 : L2;
    int curr_level_ways = (level == 1) ? L1_NUM_WAYS : L2_NUM_WAYS;
    int curr_set = (addr >> OFFSET_BITS) % (level == 1 ? L1_NUM_SETS : L2_NUM_SETS);
    unsigned curr_tag = (addr >> (OFFSET_BITS + (level == 1 ? L1_SET_BITS : L2_SET_BITS)));
    
    for (int i = 0; i < curr_level_ways; i++)
    {
        if (curr_level[i][curr_set].invalid)
        {
            curr_level[i][curr_set].addr = addr;
            curr_level[i][curr_set].tag = curr_tag;
            curr_level[i][curr_set].invalid = false;
            curr_level[i][curr_set].dirty = false;
            update_LRU(level, i, curr_set);
            return;
        }
    }
}

/**
 * @brief Writes data to the cache.
 *
 * This function handles the process of writing data to the cache. It first
 * increments the instruction count. Then, it checks if the address is present
 * in the L1 cache. If it is, it handles the L1 cache hit scenario. If the
 * address is not found in the L1 cache but is found in the L2 cache, it handles
 * the L1 miss and L2 hit scenario. If the address is not found in either the L1
 * or L2 cache, it handles the L1 miss and L2 miss scenario.
 *
 * @param addr The address to write to the cache.
 */
void Cache::Write(unsigned addr)
{
    stats.instruction_count++;
    int way_L1 = find_way_in_level(addr, 1);
    int way_L2 = find_way_in_level(addr, 2);

    if (way_L1 != MISS)  // L1 hit
    {
        Cache::Write_L1_hit(addr, way_L1);
        return;
    }
    else if (way_L2 != MISS)  // L2 hit && L1 miss
    {
        Cache::Write_L1_miss_L2_hit(addr, way_L2);
        return;
    }
    else  // L2 miss && L1 miss
    {
        Cache::Write_L1_miss_L2_miss(addr);
        return;
    }
}

/**
 * @brief Handles a write operation when there is a hit in the L1 cache.
 *
 * This function updates the cache statistics, marks the block as dirty,
 * and updates the LRU (Least Recently Used) status for the given cache line.
 *
 * @param addr The address being written to.
 * @param way_L1 The way index in the L1 cache where the hit occurred.
 */
void Cache::Write_L1_hit(unsigned addr, int way_L1)
{
    // update stats
    stats.t_access += L1Cyc;
    stats.L1_hit_count++;
    
    int set_L1 = (addr >> OFFSET_BITS) % L1_NUM_SETS;
    update_LRU(1, way_L1, set_L1);
    dirty_block(1, way_L1, set_L1);    
    return;
}

/**
 * @brief Handles the scenario where a write operation results in an L1 cache miss but an L2 cache hit.
 * 
 * This function updates the cache statistics, manages the LRU (Least Recently Used) policy, and handles
 * the dirty block status for both L1 and L2 caches. Depending on the write allocation policy, it may also
 * allocate a new block in the L1 cache.
 * 
 * @param addr The address being written to.
 * @param way_L2 The way index in the L2 cache where the block is found.
 */
void Cache::Write_L1_miss_L2_hit(unsigned addr, int way_L2)
{
    // update stats
    stats.t_access += L1Cyc;
    stats.t_access += L2Cyc;
    stats.L2_hit_count++;
    
    int set_L2 = (addr >> OFFSET_BITS) % L2_NUM_SETS;
    update_LRU(2, way_L2, set_L2); 
    dirty_block(2, way_L2, set_L2);   
    if (WrAlloc == WRITE_NO_ALOOCATE)
    {
        return;
    }
    else if (WrAlloc == WRITE_ALOOCATE)
    {
        int set_L1 = (addr >> OFFSET_BITS) % L1_NUM_SETS;
        if (level_full(1, set_L1))
        {
            int victim_way = choose_victim_way(1, set_L1);
            if (L1[victim_way][set_L1].dirty)
            {
                Block& victim_block = L1[victim_way][set_L1];
                update_lower_level(victim_block.addr); // implements update LRU
            }
            invalidate_block(1, victim_way, set_L1);
        }
        add_to_level(addr, 1); // supposed to add only to L1 and to update LRU
        dirty_block(1, find_way_in_level(addr,1), set_L1);
    }
    return;
}

/**
 * @brief Handles the case when a write operation results in a miss in both L1 and L2 caches.
 * 
 * This function updates the access time statistics to account for the cycles taken by L1, L2, 
 * and main memory accesses. Depending on the write allocation policy, it either returns immediately 
 * (WRITE_NO_ALLOCATE) or allocates a new block in the cache (WRITE_ALLOCATE).
 * 
 * @param addr The address to be written to the cache.
 * 
 * If the write allocation policy is WRITE_ALLOCATE, the function:
 * - Adds the address to the cache.
 * - Marks the corresponding blocks in both L1 and L2 caches as dirty.
 * 
 * The function calculates the set and way indices for both L1 and L2 caches to locate the blocks 
 * that need to be marked as dirty.
 */
void Cache::Write_L1_miss_L2_miss(unsigned addr)
{
    // update stats
    stats.t_access += L1Cyc;
    stats.t_access += L2Cyc;
    stats.t_access += MemCyc;

    if (WrAlloc == WRITE_NO_ALOOCATE)
    {
        return;
    }
    else if (WrAlloc == WRITE_ALOOCATE)
    {
        add_to_cache(addr);

        // dirty both instances in L1 and L2
        int set_L1 = (addr >> OFFSET_BITS) % L1_NUM_SETS;
        int way_L1 = find_way_in_level(addr, 1);
        int set_L2 = (addr >> OFFSET_BITS) % L2_NUM_SETS;
        int way_L2 = find_way_in_level(addr, 2);
        dirty_block(1, way_L1, set_L1);
        dirty_block(2, way_L2, set_L2);
        return;
    }
    return;
    
}

/**
 * @brief Calculates the L1 cache miss rate.
 *
 * This function computes the miss rate for the L1 cache by subtracting the 
 * ratio of L1 cache hits to the total number of instructions from 1.
 *
 * @return The L1 cache miss rate as a double.
 */
double Cache::get_L1_miss_rate()
{
    return 1 - (static_cast<double>(stats.L1_hit_count) / static_cast<double>(stats.instruction_count));
}

/**
 * @brief Calculate the L2 cache miss rate.
 *
 * This function computes the miss rate for the L2 cache by subtracting the 
 * ratio of L2 hits to the total number of accesses (excluding L1 hits) from 1.
 *
 * @return The L2 cache miss rate as a double.
 */
double Cache::get_L2_miss_rate()
{
    return 1 - (static_cast<double>(stats.L2_hit_count) / (static_cast<double>(stats.instruction_count) - static_cast<double>(stats.L1_hit_count)));
}

/**
 * @brief Calculates the average access time for the cache.
 *
 * This function computes the average access time by considering the miss rates
 * of both L1 and L2 caches and the respective cycle times for L1, L2, and memory.
 *
 * @return The average access time as a double.
 */
double Cache::get_avg_access_time()
{
    double l1_miss_rate = 1 - (static_cast<double>(stats.L1_hit_count) / static_cast<double>(stats.instruction_count));
    double l2_miss_rate = 1 - (static_cast<double>(stats.L2_hit_count) / (static_cast<double>(stats.instruction_count) - static_cast<double>(stats.L1_hit_count)));
    return L1Cyc + l1_miss_rate * (L2Cyc + l2_miss_rate * MemCyc);
}

/**
 * @brief Prints the current state of the cache to a debug file.
 * 
 * This function outputs the contents of both L1 and L2 caches, including
 * address, validity, dirty bit, tag, and LRU information for each way and set.
 * It also prints cache statistics such as instruction count, L1 hit count,
 * L2 hit count, and total access time.
 * 
 * @param iteration The current iteration number to be printed in the debug file.
 */
void Cache::print_cache(int iteration)
{
    std::ofstream debugFile("cache_debug_output.txt", std::ios::app); // Open file in append mode
    debugFile << "ITERATION: " << iteration << std::endl;
    debugFile << "L1:" << std::endl;
    for (int i = 0; i < L1_NUM_WAYS; i++)
    {
        debugFile << "Way " << i << ": ";
        for (int j = 0; j < L1_NUM_SETS; j++)
        {
            debugFile << "AD:_" << std::hex << L1[i][j].addr << std::dec << "   ";    // Address
            debugFile << "V:_" << !L1[i][j].invalid << "    "; // Valid
            debugFile << "D:_" << L1[i][j].dirty << "    ";    // Dirty
            debugFile << "T:_" << L1[i][j].tag << "    ";    // Tag
            debugFile << "LR:_" << L1[i][j].LRU << " | ";   // LRU
        }
        debugFile << std::endl;
    }
    debugFile << "----------------------------------------" << std::endl;
    debugFile << "L2:" << std::endl;
    for (int i = 0; i < L2_NUM_WAYS; i++)
    {
        debugFile << "Way " << i << ": ";
        for (int j = 0; j < L2_NUM_SETS; j++)
        {
            debugFile << "AD:_" << std::hex << L2[i][j].addr << std::dec << "    ";    // Address
            debugFile << "V:_" << !L2[i][j].invalid << "    "; // Valid
            debugFile << "D:_" << L2[i][j].dirty << "    ";    // Dirty
            debugFile << "T:_" << L2[i][j].tag << "    ";    // Tag
            debugFile << "LR:_" << L2[i][j].LRU << " | ";   // LRU
        }
        debugFile << std::endl;
    }
    debugFile << endl;
    debugFile << "Stats:" << std::endl;
    debugFile << "ins cnt:_" << stats.instruction_count << "    ";
    debugFile << "L1_H:_" << stats.L1_hit_count << "    ";
    debugFile << "L2_H:_" << stats.L2_hit_count << "    ";
    debugFile << "tot time:_" << stats.t_access << endl;
    debugFile << endl;
    debugFile << "########################################" << endl;
    debugFile << "########################################" << endl;
    debugFile << endl;
    debugFile << endl;
    debugFile.close(); // Close the file
}

/* 	
	========================================================================================== 
	==========================================================================================
	==========================================================================================
	==========================================================================================
*/

int main(int argc, char **argv) {

	if (argc < 19) {
		cerr << "Not enough arguments" << endl;
		return 0;
	}

	// Get input arguments

	// File
	// Assuming it is the first argument
	char* fileString = argv[1];
	ifstream file(fileString); //input file stream
	string line;
	if (!file || !file.good()) {
		// File doesn't exist or some other error
		cerr << "File not found" << endl;
		return 0;
	}

	unsigned MemCyc = 0, BSize = 0, L1Size = 0, L2Size = 0, L1Assoc = 0,
			L2Assoc = 0, L1Cyc = 0, L2Cyc = 0, WrAlloc = 0;

	for (int i = 2; i < 19; i += 2) {
		string s(argv[i]);
		if (s == "--mem-cyc") {
			MemCyc = atoi(argv[i + 1]);
		} else if (s == "--bsize") {
			BSize = atoi(argv[i + 1]);
		} else if (s == "--l1-size") {
			L1Size = atoi(argv[i + 1]);
		} else if (s == "--l2-size") {
			L2Size = atoi(argv[i + 1]);
		} else if (s == "--l1-cyc") {
			L1Cyc = atoi(argv[i + 1]);
		} else if (s == "--l2-cyc") {
			L2Cyc = atoi(argv[i + 1]);
		} else if (s == "--l1-assoc") {
			L1Assoc = atoi(argv[i + 1]);
		} else if (s == "--l2-assoc") {
			L2Assoc = atoi(argv[i + 1]);
		} else if (s == "--wr-alloc") {
			WrAlloc = atoi(argv[i + 1]);
		} else {
			cerr << "Error in arguments" << endl;
			return 0;
		}
	}

	Cache cache(MemCyc, BSize, L1Size, L2Size, L1Assoc, L2Assoc, L1Cyc, L2Cyc, WrAlloc);

    // AVIV print-debug #2. comment when running tests
    // int i = 0;
    // remove("cache_debug_output.txt");

	while (getline(file, line)) {

		stringstream ss(line);
		string address;
		char operation = 0; // read (R) or write (W)
		if (!(ss >> operation >> address)) {
			// Operation appears in an Invalid format
			cout << "Command Format error" << endl;
			return 0;
		}

		// DEBUG - remove this line
		// cout << "operation: " << operation;

		string cutAddress = address.substr(2); // Removing the "0x" part of the address

		// DEBUG - remove this line
		// cout << ", address (hex)" << cutAddress;

		unsigned long int num = 0;
		num = strtoul(cutAddress.c_str(), NULL, 16);

		if (operation == 'r') 
		{
			cache.Read(num);
		} else if (operation == 'w')
		{
			cache.Write(num);
		}
        
        // AVIV print-debug #2. comment when running tests
        // cache.print_cache(i);
        // i++;
		
        // DEBUG - remove this line
		// cout << " (dec) " << num << endl;
	}
	double L1MissRate = cache.get_L1_miss_rate();
	double L2MissRate = cache.get_L2_miss_rate();
	double avgAccTime = cache.get_avg_access_time();

	printf("L1miss=%.03f ", L1MissRate);
	printf("L2miss=%.03f ", L2MissRate);
	printf("AccTimeAvg=%.03f\n", avgAccTime);

	return 0;
}
