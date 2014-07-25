#ifndef RTREE_H
#define RTREE_H



// NOTE These next few lines may be win32 specific, you may need to modify them to compile on other platform
#define _CRT_SECURE_NO_DEPRECATE
#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <vector>

using namespace std;

#define RBSIZE 10 //we give a standar size for Reservation Buffer
#define IUNITPACKAGE 5 //number of IndexUnits in a block
#define PAGESIZE 1000 //size of every block in disk

#define ASSERT assert // RTree uses ASSERT( condition )
#ifndef Min
  #define Min __min 
#endif //Min
#ifndef Max
  #define Max __max 
#endif //Max

//
// RTree.h
//

#define RTREE_TEMPLATE template<class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES>
#define RTREE_QUAL RTree<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>

#define RTREE_DONT_USE_MEMPOOLS // This version does not contain a fixed memory allocator, fill in lines with EXAMPLE to implement one.
#define RTREE_USE_SPHERICAL_VOLUME // Better split classification, may be slower on some systems

class RTFileStream ;//: public RTree<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES> // File I/O helper class, look below for implementation and notes.

struct int_2   
{
    int i_0;
    int i_1;
};
std:: multimap <int, int_2> NodeTransTable;

int RBcount =0;		//counter which show us the next free place int ResBuffer
int page_id = 0;   //counter which show us the next free disk page 
int  UniqueId = 0;  //node unique id
int rt_page_id = 0;  
int update_flag = 0;
int IUcount = 0;

fstream mfile ;
fstream mfile2;
FILE* file;
int j;

/// class RTree
/// Implementation of RTree, a multidimensional bounding rectangle tree.
/// Example usage: For a 3-dimensional tree use RTree<Object*, float, 3> myTree;
///
/// DATATYPE Referenced data, should be int, void*, obj* etc. no larger than sizeof<void*> and simple type
/// ELEMTYPE Type of element such as int or float
/// NUMDIMS Number of dimensions such as 2 or 3  
/// ELEMTYPEREAL Type of element that allows fractional and large values such as float or double, for use in volume calcs
///
/// NOTES: Inserting and removing data requires the knowledge of its constant Minimal Bounding Rectangle.
///        This version uses new/delete for nodes, I recommend using a fixed size allocator for efficiency.
///        Instead of using a callback function for returned results, I recommend and efficient pre-sized, grow-only memory
///        array similar to MFC CArray or STL Vector for returning search query result.

///


	  
template<class DATATYPE, class ELEMTYPE, int NUMDIMS, 
         class ELEMTYPEREAL = ELEMTYPE, int TMAXNODES = 8, int TMINNODES = TMAXNODES / 2>
class RTree
{


protected: 

  struct Node;  // Fwd decl.  Used by other internal structs and iterator
  
public:

  // These constant must be declared after Branch and before Node struct
  // Stuck up here for MSVC 6 compiler.  NSVC .NET 2003 is much happier.
  enum
  {
    MAXNODES = TMAXNODES,                         ///< Max elements in node
    MINNODES = TMINNODES,                         ///< Min elements in node
  };


public:


  RTree();
  virtual ~RTree();
  

  /// Insert entry
  /// \param a_min Min of bounding rect
  /// \param a_max Max of bounding rect
  /// \param a_dataId Positive Id of data.  Maybe zero, but negative numbers not allowed.
  void Insert(const ELEMTYPE a_min[NUMDIMS], const ELEMTYPE a_max[NUMDIMS], const DATATYPE& a_dataId);
  
  /// Remove entry
  /// \param a_min Min of bounding rect
  /// \param a_max Max of bounding rect
  /// \param a_dataId Positive Id of data.  Maybe zero, but negative numbers not allowed.
  void Remove(const ELEMTYPE a_min[NUMDIMS], const ELEMTYPE a_max[NUMDIMS], const DATATYPE& a_dataId);
  
  /// Find all within search rectangle
  /// \param a_min Min of search bounding rect
  /// \param a_max Max of search bounding rect
  /// \param a_searchResult Search result array.  Caller should set grow size. Function will reset, not append to array.
  /// \param a_resultCallback Callback function to return result.  Callback should return 'true' to continue searching
  /// \param a_context User context to pass as parameter to a_resultCallback
  /// \return Returns the number of entries found
  int Search(const ELEMTYPE a_min[NUMDIMS], const ELEMTYPE a_max[NUMDIMS], bool __cdecl a_resultCallback(DATATYPE a_data, void* a_context), void* a_context);
  
  /// Remove all entries from tree
  void RemoveAll();

  /// Count the data elements in this container.  This is slow as no internal counter is maintained.
  int Count();

  /// Load tree contents from file
  bool Load(const char* a_fileName);
  /// Load tree contents from stream
  bool Load(RTFileStream& a_stream);

  
  /// Save tree contents to file
  bool Save(const char* a_fileName);
  /// Save tree contents to stream
  bool Save(RTFileStream& a_stream);
 //  void mysave( int *s_id,int min[],int max[]);//, Node *curr, Node *parent, Node *child);

  /// Iterator is not remove safe.
  class Iterator
  {
  private:
  
    enum { MAX_STACK = 32 }; //  Max stack size. Allows almost n^32 where n is number of branches in node
    
    struct StackElement
    {
      Node* m_node;
      int m_branchIndex;
    };
    
  public:
  
    Iterator()                                    { Init(); }

    ~Iterator()                                   { }
    
    /// Is iterator invalid
    bool IsNull()                                 { return (m_tos <= 0); }

    /// Is iterator pointing to valid data
    bool IsNotNull()                              { return (m_tos > 0); }

    /// Access the current data element. Caller must be sure iterator is not NULL first.
    DATATYPE& operator*()
    {
      ASSERT(IsNotNull());
      StackElement& curTos = m_stack[m_tos - 1];
      return curTos.m_node->m_branch[curTos.m_branchIndex].m_data;
    } 

    /// Access the current data element. Caller must be sure iterator is not NULL first.
    const DATATYPE& operator*() const
    {
      ASSERT(IsNotNull());
      StackElement& curTos = m_stack[m_tos - 1];
      return curTos.m_node->m_branch[curTos.m_branchIndex].m_data;
    } 

    /// Find the next data element
    bool operator++()                             { return FindNextData(); }

    /// Get the bounds for this node
    void GetBounds(ELEMTYPE a_min[NUMDIMS], ELEMTYPE a_max[NUMDIMS])
    {
      ASSERT(IsNotNull());
      StackElement& curTos = m_stack[m_tos - 1];
      Branch& curBranch = curTos.m_node->m_branch[curTos.m_branchIndex];
      
      for(int index = 0; index < NUMDIMS; ++index)
      {
        a_min[index] = curBranch.m_rect.m_min[index];
        a_max[index] = curBranch.m_rect.m_max[index];
      }
    }

  private:
  
    /// Reset iterator
    void Init()                                   { m_tos = 0; }

    /// Find the next data element in the tree (For internal use only)
    bool FindNextData()
    {
      for(;;)
      {
        if(m_tos <= 0)
        {
          return false;
        }
        StackElement curTos = Pop(); // Copy stack top cause it may change as we use it

        if(curTos.m_node->IsLeaf())
        {
          // Keep walking through data while we can
          if(curTos.m_branchIndex+1 < curTos.m_node->m_count)
          {
            // There is more data, just point to the next one
            Push(curTos.m_node, curTos.m_branchIndex + 1);
            return true;
          }
          // No more data, so it will fall back to previous level
        }
        else
        {
          if(curTos.m_branchIndex+1 < curTos.m_node->m_count)
          {
            // Push sibling on for future tree walk
            // This is the 'fall back' node when we finish with the current level
            Push(curTos.m_node, curTos.m_branchIndex + 1);
          }
          // Since cur node is not a leaf, push first of next level to get deeper into the tree
          Node* nextLevelnode = curTos.m_node->m_branch[curTos.m_branchIndex].m_child;
          Push(nextLevelnode, 0);
          
          // If we pushed on a new leaf, exit as the data is ready at TOS
          if(nextLevelnode->IsLeaf())
          {
            return true;
          }
        }
      }
    }

    /// Push node and branch onto iteration stack (For internal use only)
    void Push(Node* a_node, int a_branchIndex)
    {
      m_stack[m_tos].m_node = a_node;
      m_stack[m_tos].m_branchIndex = a_branchIndex;
      ++m_tos;
      ASSERT(m_tos <= MAX_STACK);
    }
    
    /// Pop element off iteration stack (For internal use only)
    StackElement& Pop()
    {
      ASSERT(m_tos > 0);
      --m_tos;
      return m_stack[m_tos];
    }

    StackElement m_stack[MAX_STACK];              ///< Stack as we are doing iteration instead of recursion
    int m_tos;                                    ///< Top Of Stack index
  
    friend RTree; // Allow hiding of non-public functions while allowing manipulation by logical owner
  };
   
  /// Get 'first' for iteration
  void GetFirst(Iterator& a_it)
  {
    a_it.Init();
    Node* first = m_root;
    while(first)
    {
      if(first->IsInternalNode() && first->m_count > 1)
      {
        a_it.Push(first, 1); // Descend sibling branch later
      }
      else if(first->IsLeaf())
      {
        if(first->m_count)
        {
          a_it.Push(first, 0);
        }
        break;
      }
      first = first->m_branch[0].m_child;
    }
  }  

  /// Get Next for iteration
  void GetNext(Iterator& a_it)                    { ++a_it; }

  /// Is iterator NULL, or at end?
  bool IsNull(Iterator& a_it)                     { return a_it.IsNull(); }

  /// Get object at iterator position
  DATATYPE& GetAt(Iterator& a_it)                 { return *a_it; }

protected:

  /// Minimal bounding rectangle (n-dimensional)
  struct Rect
  {
    ELEMTYPE m_min[NUMDIMS];                      ///< Min dimensions of bounding box 
    ELEMTYPE m_max[NUMDIMS];                      ///< Max dimensions of bounding box 
  };

  /// May be data or may be another subtree
  /// The parents level determines this.
  /// If the parents level is 0, then this is data
  struct Branch
  {
	int iu_count;								  //Node place inside iunit
    Rect m_rect;                                  ///< Bounds
    union
    {
      Node* m_child;                              ///< Child node
      DATATYPE m_data;                            ///< Data Id or Ptr
	  
    };
  };

  /// Node for each branch level
  struct Node
  {
    bool IsInternalNode()                         { return (m_level > 0); } // Not a leaf, but a internal node
     bool IsLeaf()                                 { return (m_level == 0); } // A leaf, contains data
    
    int unique_id;								  //Node id
    int m_count;                                  ///< Count
    int m_level;                                  ///< Leaf is zero, others positive
    Branch m_branch[MAXNODES];                    ///< Branch
  }; 
  
  /// A link list of nodes for reinsertion after a delete operation
  struct ListNode
  {
    ListNode* m_next;                             ///< Next in list
    Node* m_node;                                 ///< Node
  };

  
   struct IndexUnit{
	 
	 struct Rect minimal_bounding_box;
	 DATATYPE *id;
	 char *op_flag;
	 Node *child_ptr;
	 Node *parent_ptr;
	 DATATYPE *data_ptr;
	 Node *IUNode;
	 int iu_id;

	 void FillIndexUnit();

	 }IUnit;

  std::list< IndexUnit> mylist;

  public:
	 
  //Reservatinon Buffer
  struct ResBuffer
  {
	ELEMTYPE min[NUMDIMS];
	ELEMTYPE max[NUMDIMS];
	//struct Rect *RBrect;

	DATATYPE *id;
	char *op;
	int count;
	void FillResBuffer();
  }RB[RBSIZE];



 //function which fill the Reservation Buffer  
  void FillResBuffer(Rect *_rect, DATATYPE c, char d, int count){

	 for(int i=0;i<NUMDIMS;i++){
		 RB[count].min[i] =  _rect->m_min[i];
		 RB[count].max[i] =  _rect->m_max[i];

	 }
	
	 RB[count].id = (int*)c;
	 RB[count].op = (char*)d;
 
  }

  //When ResBuffer is full then Reset it.
  void ResetResBuffer(void)
  {

		for(int i=0;i<RBSIZE;i++){
		
			for(int j=0;j<NUMDIMS;j++){
				RB[i].min[j] = NULL;
				RB[i].max[j] = NULL;
			}
			RB[i].id = NULL;
			RB[i].op = NULL;
		}
	// free(RB);
	 RBcount = 0;

	}

  Node *nodebuffer[2];

 

  void FillIndexUnit(ResBuffer RB[],int RBindex)
  {

	IUnit.id = RB[RBindex].id;
	IUnit.op_flag = RB[RBindex].op;


	for(int j=0;j<NUMDIMS;j++){

		IUnit.minimal_bounding_box.m_min[j]= RB[RBindex].min[j];
	
		IUnit.minimal_bounding_box.m_max[j]= RB[RBindex].max[j];

	}
	IUnit.IUNode = nodebuffer[1];
	IUnit.parent_ptr = nodebuffer[0];
	IUnit.child_ptr = NULL;
	IUnit.iu_id = (int)IUnit.IUNode->unique_id;//}
	mylist.push_back(IUnit);
	
	if((RBindex+1)%IUNITPACKAGE == 0 && RBindex!=0){

		Save(mylist);
		mylist.clear();
		IUcount=0;
	
	 }
	
   }

  //save an index unit structure to disk and make the insertion to NodeTransTable
  void  Save(std::list<IndexUnit> &mylist)
  {
	
	mfile.open("savefile",ios::out|ios::in|std::ios::binary);

	char buffer[PAGESIZE];
	int IUtotalsize = 6*sizeof(ELEMTYPE)+sizeof(char); 
	memset(buffer,'0',PAGESIZE);
	int i = 0;
	int_2 table ;
	int flag = 0;
	int count =0;
	int list_count = 1;
	int *temp;
	

	for(std::list<IndexUnit>::iterator list_it=mylist.begin(); list_it!=mylist.end(); ++list_it)
	{
	 
	   temp = (int*)list_it->iu_id;

	   memcpy(buffer+(i*IUtotalsize),&list_it->id,sizeof(int));
	   memcpy(buffer+(i*IUtotalsize)+sizeof(ELEMTYPE),&list_it->minimal_bounding_box.m_min[0],sizeof(ELEMTYPE));
	   memcpy(buffer+(i*IUtotalsize)+2*sizeof(ELEMTYPE),&list_it->minimal_bounding_box.m_min[1],sizeof(ELEMTYPE));
	   memcpy(buffer+(i*IUtotalsize)+3*sizeof(ELEMTYPE),&list_it->minimal_bounding_box.m_max[0],sizeof(ELEMTYPE)); 
	   memcpy(buffer+(i*IUtotalsize)+4*sizeof(ELEMTYPE),&list_it->minimal_bounding_box.m_max[1],sizeof(ELEMTYPE));
	   memcpy(buffer+(i*IUtotalsize)+5*sizeof(ELEMTYPE),&list_it->op_flag,sizeof(char));
	
	   list_count++;
	
	   memcpy(buffer+(i*IUtotalsize)+5*sizeof(ELEMTYPE)+sizeof(char),&temp,sizeof(int));
	
	   table.i_0 = page_id;
	   table.i_1 = i;
	   NodeTransTable.insert(std::pair<int,int_2>((const int) temp,table));

		++i;
	 }
  
	 update_flag = 0;
	 int seek = page_id*PAGESIZE;

	 ++page_id;

     mfile.seekp(seek,std::ios_base::beg);
     mfile.write (buffer, PAGESIZE);
	
	 if (mfile.fail()){
		
		mfile.clear();
		mfile.write (buffer, PAGESIZE);
	 }

     mfile.flush();

     if(mfile)
     {
	   mfile.close();
     }

  }


  /// Variables for finding a split partition
  struct PartitionVars
  {
    int m_partition[MAXNODES+1];
    int m_total;
    int m_minFill;
    int m_taken[MAXNODES+1];
    int m_count[2];
    Rect m_cover[2];
    ELEMTYPEREAL m_area[2];

    Branch m_branchBuf[MAXNODES+1];
    int m_branchCount;
    Rect m_coverSplit;
    ELEMTYPEREAL m_coverSplitArea;
  }; 
 
  Node* AllocNode();
  void FreeNode(Node* a_node);
  void InitNode(Node* a_node);
  void InitRect(Rect* a_rect);
  bool InsertRectRec(Rect* a_rect, const DATATYPE& a_id, Node* a_node, Node** a_newNode, int a_level, int iucount, int interval);
  bool InsertRect(Rect* a_rect, const DATATYPE& a_id, Node** a_root, int a_level, int iucount, int interval);
  Rect NodeCover(Node* a_node);
  bool AddBranch(Branch* a_branch, Node* a_node, Node** a_newNode);
  void DisconnectBranch(Node* a_node, int a_index);
  int PickBranch(Rect* a_rect, Node* a_node);
  Rect CombineRect(Rect* a_rectA, Rect* a_rectB);
  void SplitNode(Node* a_node, Branch* a_branch, Node** a_newNode);
  ELEMTYPEREAL RectSphericalVolume(Rect* a_rect);
  ELEMTYPEREAL RectVolume(Rect* a_rect);
  ELEMTYPEREAL CalcRectVolume(Rect* a_rect);
  void GetBranches(Node* a_node, Branch* a_branch, PartitionVars* a_parVars);
  void ChoosePartition(PartitionVars* a_parVars, int a_minFill);
  void LoadNodes(Node* a_nodeA, Node* a_nodeB, PartitionVars* a_parVars);
  void InitParVars(PartitionVars* a_parVars, int a_maxRects, int a_minFill);
  void PickSeeds(PartitionVars* a_parVars);
  void Classify(int a_index, int a_group, PartitionVars* a_parVars);
  bool RemoveRect(Rect* a_rect, const DATATYPE& a_id, Node** a_root);
  bool RemoveRectRec(Rect* a_rect, const DATATYPE& a_id, Node* a_node, ListNode** a_listNode);
  ListNode* AllocListNode();
  void FreeListNode(ListNode* a_listNode);
  bool Overlap(Rect* a_rectA, Rect* a_rectB);
  void ReInsert(Node* a_node, ListNode** a_listNode);
  bool Search(Node* a_node, Rect* a_rect, int& a_foundCount, bool __cdecl a_resultCallback(DATATYPE a_data, void* a_context), void* a_context);
  void RemoveAllRec(Node* a_node);
  void Reset();
  void CountRec(Node* a_node, int& a_count);
  bool SaveRec(Node* a_node, RTFileStream& a_stream);
  bool LoadRec(Node* a_node, RTFileStream& a_stream);
  bool SaveRecTwo(Node* a_node, RTFileStream& a_stream);
  bool SaveRecThree(Node* a_node, RTFileStream& a_stream);
  bool LoadRecTwo(Node* a_node, RTFileStream& a_stream);
  bool LoadRecThree(Node* a_node, RTFileStream& a_stream);
  void ConstractNode(Node** a_node, int node_id);
  Node* m_root;                                    ///< Root of tree
  ELEMTYPEREAL m_unitSphereVolume;                 ///< Unit sphere constant for required number of dimensions
  

};


// Because there is not stream support, this is a quick and dirty file I/O helper.
// Users will likely replace its usage with a Stream implementation from their favorite API.
//template<class DATATYPE, class ELEMTYPE, int NUMDIMS, 
//         class ELEMTYPEREAL = ELEMTYPE, int TMAXNODES = 8, int TMINNODES = TMAXNODES / 2>
class RTFileStream //: public RTree<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>
{

  FILE* m_file;
  
public:
	
  RTFileStream()
  {
    m_file = NULL;
  }

  ~RTFileStream()
  {
    Close();
  }

  bool OpenRead(const char* a_fileName)
  {
    m_file = fopen(a_fileName, "rb");
    if(!m_file)
    {
      return false;
    }
    return true;
  }

  bool OpenWrite(const char* a_fileName)
  {
    m_file = fopen(a_fileName, "wb");
    if(!m_file)
    {
      return false;
    }
    return true;
  }

  void Close()
  {
    if(m_file)
    {
      fclose(m_file);
      m_file = NULL;
    }
  }

  template< typename TYPE >
  size_t Write(const TYPE& a_value)
  {
    ASSERT(m_file);
    return fwrite((void*)&a_value, sizeof(a_value), 1, m_file);
  }

  template< typename TYPE >
  size_t WriteArray(const TYPE* a_array, int a_count)
  {
    ASSERT(m_file);
    return fwrite((void*)a_array, sizeof(TYPE) * a_count, 1, m_file);
  }

  template< typename TYPE >
  size_t Read(TYPE& a_value)
  {
    ASSERT(m_file);
    return fread((void*)&a_value, sizeof(a_value), 1, m_file);
  }

  template< typename TYPE >
  size_t ReadArray(TYPE* a_array, int a_count)
  {
    ASSERT(m_file);
    return fread((void*)a_array, sizeof(TYPE) * a_count, 1, m_file);
  }
 
};


RTREE_TEMPLATE
RTREE_QUAL::RTree()
{
  ASSERT(MAXNODES > MINNODES);
  ASSERT(MINNODES > 0);


  // We only support machine word size simple data type eg. integer index or object pointer.
  // Since we are storing as union with non data branch
  ASSERT(sizeof(DATATYPE) == sizeof(void*) || sizeof(DATATYPE) == sizeof(int));

  // Precomputed volumes of the unit spheres for the first few dimensions
  const float UNIT_SPHERE_VOLUMES[] = {
    0.000000f, 2.000000f, 3.141593f, // Dimension  0,1,2
    4.188790f, 4.934802f, 5.263789f, // Dimension  3,4,5
    5.167713f, 4.724766f, 4.058712f, // Dimension  6,7,8
    3.298509f, 2.550164f, 1.884104f, // Dimension  9,10,11
    1.335263f, 0.910629f, 0.599265f, // Dimension  12,13,14
    0.381443f, 0.235331f, 0.140981f, // Dimension  15,16,17
    0.082146f, 0.046622f, 0.025807f, // Dimension  18,19,20 
  };

  m_root = AllocNode();
  m_root->m_level = 0;

  m_unitSphereVolume = (ELEMTYPEREAL)UNIT_SPHERE_VOLUMES[NUMDIMS];
  
}


RTREE_TEMPLATE
RTREE_QUAL::~RTree()
{
  Reset(); // Free, or reset node memory
}


RTREE_TEMPLATE
void RTREE_QUAL::Insert(const ELEMTYPE a_min[NUMDIMS], const ELEMTYPE a_max[NUMDIMS], const DATATYPE& a_dataId)
{
	

#ifdef _DEBUG
  for(int index=0; index<NUMDIMS; ++index)
  {
    ASSERT(a_min[index] <= a_max[index]);
  }
#endif //_DEBUG

  Rect rect;
 
  for(int axis=0; axis<NUMDIMS; ++axis)
  { 
    rect.m_min[axis] = a_min[axis];
    rect.m_max[axis] = a_max[axis];
  }
 
  if (RBcount<RBSIZE)
  {
	  FillResBuffer(&rect,a_dataId,'i',RBcount);
  }

  else
  {
	  ResetResBuffer();
	  }

  RBcount++;

  if(RBcount==RBSIZE)
  {
	IUcount = 0;
		
	for(int i=0;i<RBSIZE;i++){

		for(int j=0; j<NUMDIMS; ++j)
		{
			rect.m_min[j] = RB[i].min[j];
			rect.m_max[j] = RB[i].max[j];
		 } 
		const int &Id =  (DATATYPE&)RB[i].id;
		
		nodebuffer[0]=NULL;
		nodebuffer[1]=NULL;

		if(RB[i].op ==  (char*)'i')
		{		
		  InsertRect(&rect, Id, &m_root, 0, 0, 0);
		  FillIndexUnit(RB,i);
		}
		else if( RB[i].op ==  (char*)'r')
		{
		  RemoveRect(&rect, Id, &m_root);
		  FillIndexUnit(RB,i);
		}
	 }

	RBcount=0;

	}
}


RTREE_TEMPLATE
void RTREE_QUAL::Remove(const ELEMTYPE a_min[NUMDIMS], const ELEMTYPE a_max[NUMDIMS], const DATATYPE& a_dataId)
{
  #ifdef _DEBUG
    for(int index=0; index<NUMDIMS; ++index)
	{
    ASSERT(a_min[index] <= a_max[index]);
	}
  #endif //_DEBUG

  Rect rect;
 
  for(int axis=0; axis<NUMDIMS; ++axis)
  {
    rect.m_min[axis] = a_min[axis];
    rect.m_max[axis] = a_max[axis];
  }

  if (RBcount<RBSIZE)
  {
	  FillResBuffer(&rect,a_dataId,'r',RBcount);
  }

  else
  {
	  ResetResBuffer();
  }

  RBcount++;

  if(RBcount==RBSIZE)
  {
	IUcount = 0;

	for(int i=0;i<RBSIZE;i++)
	{
		 for(int j=0; j<NUMDIMS; ++j)
		 {
			rect.m_min[j] = RB[i].min[j];
			rect.m_max[j] = RB[i].max[j];
		 }
		 
		const int &Id =  (DATATYPE&)RB[i].id;
		nodebuffer[0]=NULL;
		nodebuffer[1]=NULL;

		if(RB[i].op == (char*)'i')
		{
			
			InsertRect(&rect, Id, &m_root, 0, 0, 0);
			FillIndexUnit(RB,i);
		}

		else if(RB[i].op == (char*)'r')
		{
			RemoveRect(&rect, Id, &m_root);
			FillIndexUnit(RB,i);
		}
		

	}

	RBcount=0;
  }
}


RTREE_TEMPLATE
int RTREE_QUAL::Search(const ELEMTYPE a_min[NUMDIMS], const ELEMTYPE a_max[NUMDIMS], bool __cdecl a_resultCallback(DATATYPE a_data, void* a_context), void* a_context)
{
#ifdef _DEBUG
  for(int index=0; index<NUMDIMS; ++index)
  {
    ASSERT(a_min[index] <= a_max[index]);
  }
#endif //_DEBUG

  Rect rect;
  
  for(int axis=0; axis<NUMDIMS; ++axis)
  {
    rect.m_min[axis] = a_min[axis];
    rect.m_max[axis] = a_max[axis];
  }

  // NOTE: May want to return search result another way, perhaps returning the number of found elements here.

  int foundCount = 0;
  Search(m_root, &rect, foundCount, a_resultCallback, a_context);

  return foundCount;
}


RTREE_TEMPLATE
int RTREE_QUAL::Count()
{
  int count = 0;
  CountRec(m_root, count);
  
  return count;
}



RTREE_TEMPLATE
void RTREE_QUAL::CountRec(Node* a_node, int& a_count)
{
  if(a_node->IsInternalNode())  // not a leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      CountRec(a_node->m_branch[index].m_child, a_count);
    }
  }
  else // A leaf node
  {
    a_count += a_node->m_count;
  }
}


RTREE_TEMPLATE
bool RTREE_QUAL::Load(const char* a_fileName)
{
  RemoveAll(); // Clear existing tree

  RTFileStream stream;
  if(!stream.OpenRead(a_fileName))
  {
    return false;
  }

  bool result = Load(stream);
  
  stream.Close();

  return result;
};



RTREE_TEMPLATE
bool RTREE_QUAL::Load(RTFileStream& a_stream)
{
  // Write some kind of header
  int _dataFileId = ('R'<<0)|('T'<<8)|('R'<<16)|('E'<<24);
  int _dataSize = sizeof(DATATYPE);
  int _dataNumDims = NUMDIMS;
  int _dataElemSize = sizeof(ELEMTYPE);
  int _dataElemRealSize = sizeof(ELEMTYPEREAL);
  int _dataMaxNodes = TMAXNODES;
  int _dataMinNodes = TMINNODES;

  int dataFileId = 0;
  int dataSize = 0;
  int dataNumDims = 0;
  int dataElemSize = 0;
  int dataElemRealSize = 0;
  int dataMaxNodes = 0;
  int dataMinNodes = 0;

  a_stream.Read(dataFileId);
  a_stream.Read(dataSize);
  a_stream.Read(dataNumDims);
  a_stream.Read(dataElemSize);
  a_stream.Read(dataElemRealSize);
  a_stream.Read(dataMaxNodes);
  a_stream.Read(dataMinNodes);

  bool result = false;

  // Test if header was valid and compatible
  if(    (dataFileId == _dataFileId) 
      && (dataSize == _dataSize) 
      && (dataNumDims == _dataNumDims) 
      && (dataElemSize == _dataElemSize) 
      && (dataElemRealSize == _dataElemRealSize) 
      && (dataMaxNodes == _dataMaxNodes) 
      && (dataMinNodes == _dataMinNodes) 
    )
  {
    // Recursively load tree
    result = LoadRecThree(m_root, a_stream);
  }

  return result;
}

RTREE_TEMPLATE
bool RTREE_QUAL::LoadRec(Node* a_node, RTFileStream& a_stream)
{


  a_stream.Read(a_node->m_level);
  a_stream.Read(a_node->m_count);

  if(a_node->IsInternalNode())  // not a leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      Branch* curBranch = &a_node->m_branch[index];

      a_stream.ReadArray(curBranch->m_rect.m_min, NUMDIMS);
      a_stream.ReadArray(curBranch->m_rect.m_max, NUMDIMS);

      curBranch->m_child = AllocNode();
      LoadRec(curBranch->m_child, a_stream);
    }
  }
  else // A leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      Branch* curBranch = &a_node->m_branch[index];

      a_stream.ReadArray(curBranch->m_rect.m_min, NUMDIMS);
      a_stream.ReadArray(curBranch->m_rect.m_max, NUMDIMS);

      a_stream.Read(curBranch->m_data);
    }
  }

  return true; // Should do more error checking on I/O operations
}

RTREE_TEMPLATE
bool RTREE_QUAL::LoadRecTwo(Node* a_node, RTFileStream& a_stream)
{

  mfile2.open("savefile2", ios::in|ios:: out);
  char buffer[PAGESIZE];

  a_stream.Read(a_node->m_level);
  a_stream.Read(a_node->m_count);

  if(a_node->IsInternalNode())  // not a leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      Branch* curBranch = &a_node->m_branch[index];

      a_stream.ReadArray(curBranch->m_rect.m_min, NUMDIMS);
      a_stream.ReadArray(curBranch->m_rect.m_max, NUMDIMS);

      curBranch->m_child = AllocNode();
      LoadRec(curBranch->m_child, a_stream);
    }
  }
  else // A leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      Branch* curBranch = &a_node->m_branch[index];

	  memset(buffer,'0',PAGESIZE);
	
	  int seek = rt_page_id*PAGESIZE;
	  mfile2.read ((char*)buffer, sizeof(buffer));
	  mfile2.seekg(seek,std::ios_base::beg); 
	  
	  memcpy(curBranch->m_rect.m_min,buffer,2*sizeof(ELEMTYPE));
	  memcpy(curBranch->m_rect.m_max,buffer+2*sizeof(int),2*sizeof(ELEMTYPE));
	  memcpy(&curBranch->m_data,buffer+4*sizeof(ELEMTYPE),sizeof(int));
	  rt_page_id++;
	
	}

	//printf("load_enter");
	if(mfile2)
   {
		mfile2.close();
    }
  }
  return true; // Should do more error checking on I/O operations
}

RTREE_TEMPLATE
bool RTREE_QUAL::LoadRecThree(Node* a_node, RTFileStream& a_stream)
{

  mfile2.open("savefile2", ios::in|ios:: out);

  a_stream.Read(a_node->m_level);  
  a_stream.Read(a_node->m_count);
  a_stream.Read(a_node->unique_id);
  
  if(a_node->IsInternalNode())  // not a leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      Branch* curBranch = &a_node->m_branch[index];

      a_stream.ReadArray(curBranch->m_rect.m_min, NUMDIMS);
      a_stream.ReadArray(curBranch->m_rect.m_max, NUMDIMS);

      curBranch->m_child = AllocNode();
      LoadRecThree(curBranch->m_child, a_stream);
    }
  }
  else // A leaf node
  {
	Node* tmp_node;
	ConstractNode(&tmp_node, a_node->unique_id);
	int index=0;

	for(int i=0;i<tmp_node->m_count;i++)
	{
		a_node->m_branch[index]=tmp_node->m_branch[index];
		index++;
	}
	
  }

  return true; // Should do more error checking on I/O operations

}

RTREE_TEMPLATE
bool RTREE_QUAL::Save(const char* a_fileName)
{
  RTFileStream stream;
  if(!stream.OpenWrite(a_fileName))
  {
    return false;
  }

  bool result = Save(stream);
  stream.Close();

  return result;
}


RTREE_TEMPLATE
bool RTREE_QUAL::Save(RTFileStream& a_stream)
{
  // Write some kind of header
  int dataFileId = ('R'<<0)|('T'<<8)|('R'<<16)|('E'<<24);
  int dataSize = sizeof(DATATYPE);
  int dataNumDims = NUMDIMS;
  int dataElemSize = sizeof(ELEMTYPE);
  int dataElemRealSize = sizeof(ELEMTYPEREAL);
  int dataMaxNodes = TMAXNODES;
  int dataMinNodes = TMINNODES;

  a_stream.Write(dataFileId);
  a_stream.Write(dataSize);
  a_stream.Write(dataNumDims);
  a_stream.Write(dataElemSize);
  a_stream.Write(dataElemRealSize);
  a_stream.Write(dataMaxNodes);
  a_stream.Write(dataMinNodes);

  // Recursively save tree
  bool result = SaveRecThree(m_root, a_stream);
  rt_page_id = 0;
  return result;
}


RTREE_TEMPLATE
bool RTREE_QUAL::SaveRec(Node* a_node, RTFileStream& a_stream)
{
  a_stream.Write(a_node->m_level);
  a_stream.Write(a_node->m_count);

  if(a_node->IsInternalNode())  // not a leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      Branch* curBranch = &a_node->m_branch[index];

      a_stream.WriteArray(curBranch->m_rect.m_min, NUMDIMS);
      a_stream.WriteArray(curBranch->m_rect.m_max, NUMDIMS);

      SaveRec(curBranch->m_child, a_stream);
    }
  }
  else // A leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      Branch* curBranch = &a_node->m_branch[index];

      a_stream.WriteArray(curBranch->m_rect.m_min, NUMDIMS);
      a_stream.WriteArray(curBranch->m_rect.m_max, NUMDIMS);
      a_stream.Write(curBranch->m_data);
    }
  }

  return true; // Should do more error checking on I/O operations
}
RTREE_TEMPLATE
bool RTREE_QUAL::SaveRecTwo(Node* a_node, RTFileStream& a_stream)
{
  mfile2.open("savefile2",ios:: out);

  char buffer[PAGESIZE];
  memset(buffer,'0',PAGESIZE);
  
  a_stream.Write(a_node->m_level);
  a_stream.Write(a_node->m_count);
  a_stream.Write(a_node->unique_id);


  if(a_node->IsInternalNode())  // not a leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      Branch* curBranch = &a_node->m_branch[index];

      a_stream.WriteArray(curBranch->m_rect.m_min, NUMDIMS);
      a_stream.WriteArray(curBranch->m_rect.m_max, NUMDIMS);

      SaveRecTwo(curBranch->m_child, a_stream);                  //pi8ano la8os
    }
  }
  else // A leaf node
  {
	 
    for(int index = 0; index < a_node->m_count; ++index)
    {

      Branch* curBranch = &a_node->m_branch[index];
	  
	  memcpy(buffer,curBranch->m_rect.m_min,2*sizeof(int));
	  memcpy(buffer+2*sizeof(int),curBranch->m_rect.m_max,2*sizeof(int));
      memcpy(buffer+4*sizeof(int),&curBranch->m_data,sizeof(int));
	  int seek = rt_page_id*PAGESIZE;
      rt_page_id++;

	  mfile2.seekp(seek,std::ios_base::beg);
      mfile2.write (buffer, PAGESIZE);
	 
      if (mfile2.fail())
	  {
		
		mfile2.clear();
		mfile2.write (buffer, PAGESIZE);
	  }

      mfile2.flush();
   }

   if(mfile2)
   {
	  mfile2.close();
    }
 }

  return true; // Should do more error checking on I/O operations
}

RTREE_TEMPLATE
bool RTREE_QUAL::SaveRecThree(Node* a_node, RTFileStream& a_stream)
{
  mfile2.open("savefile2",ios:: out);
  char buffer[PAGESIZE];
 
  memset(buffer,'0',PAGESIZE);
  
  a_stream.Write(a_node->m_level);
  a_stream.Write(a_node->m_count);
  a_stream.Write(a_node->unique_id);

  if(a_node->IsInternalNode())  // not a leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      Branch* curBranch = &a_node->m_branch[index];

      a_stream.WriteArray(curBranch->m_rect.m_min, NUMDIMS);
      a_stream.WriteArray(curBranch->m_rect.m_max, NUMDIMS);

      SaveRecThree(curBranch->m_child, a_stream);                 
    }
  }
 // else // A leaf node do nothing
  

  return true; // Should do more error checking on I/O operations
}


RTREE_TEMPLATE
void RTREE_QUAL::RemoveAll()
{
  // Delete all existing nodes
  Reset();

  m_root = AllocNode();
  m_root->m_level = 0;
 
}


RTREE_TEMPLATE
void RTREE_QUAL::Reset()
{
#ifdef RTREE_DONT_USE_MEMPOOLS
  // Delete all existing nodes
  RemoveAllRec(m_root);
  
#else // RTREE_DONT_USE_MEMPOOLS
  // Just reset memory pools.  We are not using complex types
  // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_QUAL::RemoveAllRec(Node* a_node)
{
  ASSERT(a_node);
  ASSERT(a_node->m_level >= 0);

  if(a_node->IsInternalNode()) // This is an internal node in the tree
  {
    for(int index=0; index < a_node->m_count; ++index)
    {
      RemoveAllRec(a_node->m_branch[index].m_child);
    }
  }

  FreeNode(a_node); 
 
}


RTREE_TEMPLATE
typename RTREE_QUAL::Node* RTREE_QUAL::AllocNode()
{
  Node* newNode;
#ifdef RTREE_DONT_USE_MEMPOOLS
  newNode = new Node[2];
#else // RTREE_DONT_USE_MEMPOOLS
  // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
  InitNode(newNode);
  return newNode;
}


RTREE_TEMPLATE
void RTREE_QUAL::FreeNode(Node* a_node)
{
  ASSERT(a_node);
 
#ifdef RTREE_DONT_USE_MEMPOOLS
  delete a_node;
#else // RTREE_DONT_USE_MEMPOOLS
  // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


// Allocate space for a node in the list used in DeletRect to
// store Nodes that are too empty.
RTREE_TEMPLATE
typename RTREE_QUAL::ListNode* RTREE_QUAL::AllocListNode()
{
#ifdef RTREE_DONT_USE_MEMPOOLS
  return new ListNode;
#else // RTREE_DONT_USE_MEMPOOLS
  // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_QUAL::FreeListNode(ListNode* a_listNode)
{
#ifdef RTREE_DONT_USE_MEMPOOLS
  delete a_listNode;
#else // RTREE_DONT_USE_MEMPOOLS
  // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_QUAL::InitNode(Node* a_node)
{
  a_node->m_count = 0;
  a_node->m_level = -1;
}


RTREE_TEMPLATE
void RTREE_QUAL::InitRect(Rect* a_rect)
{
  for(int index = 0; index < NUMDIMS; ++index)
  {
    a_rect->m_min[index] = (ELEMTYPE)0;
    a_rect->m_max[index] = (ELEMTYPE)0;
  }
}


// Inserts a new data rectangle into the index structure.
// Recursively descends tree, propagates splits back up.
// Returns 0 if node was not split.  Old node updated.
// If node was split, returns 1 and sets the pointer pointed to by
// new_node to point to the new node.  Old node updated to become one of two.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
// If is an internal insertion then make the changes in NodeTransTable and IndexUnit list
RTREE_TEMPLATE
bool RTREE_QUAL::InsertRectRec(Rect* a_rect, const DATATYPE& a_id, Node* a_node, Node** a_newNode, int a_level, int iucount, int interval)
{
  ASSERT(a_rect && a_node && a_newNode);
  ASSERT(a_level >= 0 && a_level <= a_node->m_level);

  int index = 0;
  Branch branch;
  Node* otherNode;
  
  // Still above level for insertion, go down tree recursively
  if(a_node->m_level > a_level)
  {
    index = PickBranch(a_rect, a_node);

    if (!InsertRectRec(a_rect, a_id, a_node->m_branch[index].m_child, &otherNode, a_level, iucount, interval))
    {
      // Child was not split
      a_node->m_branch[index].m_rect = CombineRect(a_rect, &(a_node->m_branch[index].m_rect));
	
      return false;
    }

    else // Child was split
    {
      a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
      branch.m_child = otherNode;
      branch.m_rect = NodeCover(otherNode);
	
	  bool res = AddBranch(&branch, a_node, a_newNode);

      return res;

    }
  }
  else if(a_node->m_level == a_level) // Have reached level for insertion. Add rect, split if necessary
  {
    branch.m_rect = *a_rect;
    branch.m_child = (Node*) a_id;
	
	if(a_node->unique_id < 0)
	  {
		a_node->unique_id = ++UniqueId;    
	  }
	
	if(interval==0)
	{
		branch.iu_count = ++IUcount;
		nodebuffer[1] = a_node; 
	}

	else if(interval > 0)
	{

		branch.iu_count = iucount;

		char buffer[PAGESIZE];
		mfile.open("savefile", ios::in|ios:: out|std::ios::binary); 
		int IUtotalsize = 5*sizeof(ELEMTYPE)+sizeof(int)+sizeof(char);
		int flag = 0;
		memset(buffer,'0',PAGESIZE);
		int_2 table ;
		pair<multimap <int,int_2>::iterator, multimap <int,int_2>::iterator > map_it ;
		map_it = NodeTransTable.equal_range(interval);

		for (std::multimap <int, int_2>::iterator it2 = map_it.first; it2 != map_it.second;)
		{
			int seek = (*it2).second.i_0*PAGESIZE;
			mfile.seekg(seek,std::ios_base::beg); 
			mfile.read ((char*)buffer, sizeof(buffer));	 

			ELEMTYPE load_temp[4];
			int	load_id;
			int newid=-1;
	
			memcpy(&load_id,buffer+((*it2).second.i_1* IUtotalsize),sizeof(int));
			memcpy(&load_temp[0],buffer+((*it2).second.i_1* IUtotalsize)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&load_temp[1],buffer+((*it2).second.i_1* IUtotalsize)+sizeof(ELEMTYPE)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&load_temp[2],buffer+((*it2).second.i_1* IUtotalsize)+2*sizeof(ELEMTYPE)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&load_temp[3],buffer+((*it2).second.i_1* IUtotalsize)+3*sizeof(ELEMTYPE)+sizeof(int),sizeof(ELEMTYPE));

			if((load_temp[0] == branch.m_rect.m_min[0])&&(load_temp[1] == branch.m_rect.m_min[1])
				&&(load_temp[2] ==  branch.m_rect.m_max[0])&&(load_temp[3] ==  branch.m_rect.m_max[1])&&(load_id == branch.m_data))
			{
				table.i_0 = (*it2).second.i_0;
				table.i_1 = (*it2).second.i_1; 
				std::multimap <int, int_2>::iterator itcopy = it2;
				itcopy++;
					
				NodeTransTable.erase(it2);
				it2=itcopy;
				NodeTransTable.insert(std::pair<int,int_2>((a_node)->unique_id,table));
				flag=1;
			
				break;                                                     
			}

			else ++it2;		
		}

		
		if(flag==0)
		{ 
			
				int list_count = 1;
				std::list<int_2>::iterator s_it;

				for(std::list<IndexUnit>::iterator list_it=mylist.begin(); list_it!=mylist.end(); ++list_it)
				{
					if(list_count == branch.iu_count)
					{
						list_it->iu_id = (int)(*a_newNode)->unique_id;
						flag=2;
						break;
					
					}
					list_count++;
				}
			
		}

		if((flag == 0)&&(interval==0))
		{
			nodebuffer[1] = (*a_newNode);
		}


		if(mfile)
		{
		mfile.close();
		}

	}

	bool res_else = AddBranch(&branch, a_node, a_newNode);

    return res_else;

 }

  else
  {
    // Should never occur
    ASSERT(0);
    return false;
  }
}

// Insert a data rectangle into an index structure.
// InsertRect provides for splitting the root;
// returns 1 if root was split, 0 if it was not.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
// InsertRect2 does the recursion.
RTREE_TEMPLATE
bool RTREE_QUAL:: InsertRect(Rect* a_rect, const DATATYPE& a_id, Node** a_root, int a_level, int iucount, int interval)
{
  ASSERT(a_rect && a_root);
  ASSERT(a_level >= 0 && a_level <= (*a_root)->m_level);
 #ifdef _DEBUG
  for(int index=0; index < NUMDIMS; ++index)
  {
    ASSERT(a_rect->m_min[index] <= a_rect->m_max[index]);
  }
 #endif //_DEBUG  

  Node* newRoot;
  Node* newNode;
  Branch branch;

	if(InsertRectRec(a_rect, a_id, *a_root, &newNode, a_level, iucount, interval))  // Root split
	{
	 newRoot = AllocNode();  // Grow tree taller and new root
	 newRoot->m_level = (*a_root)->m_level + 1;
	 branch.m_rect = NodeCover(*a_root);
	 branch.m_child = *a_root;
	 AddBranch(&branch, newRoot, NULL);
	 branch.m_rect = NodeCover(newNode);
     branch.m_child = newNode;
     AddBranch(&branch, newRoot, NULL);
     *a_root = newRoot;

    return true;
	}

  return false;
}


// Find the smallest rectangle that includes all rectangles in branches of a node.
RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::NodeCover(Node* a_node)
{
  ASSERT(a_node);
 
  int firstTime = true;
  Rect rect;
  InitRect(&rect);
 
  for(int index = 0; index < a_node->m_count; ++index)
  {
    if(firstTime)
    {
      rect = a_node->m_branch[index].m_rect;
      firstTime = false;
    }
    else
    {
      rect = CombineRect(&rect, &(a_node->m_branch[index].m_rect));
    }
  }
  
  return rect;
}


// Add a branch to a node.  Split the node if necessary.
// Returns 0 if node not split.  Old node updated.
// Returns 1 if node split, sets *new_node to address of new node.
// Old node updated, becomes one of two.
RTREE_TEMPLATE
bool RTREE_QUAL::AddBranch(Branch* a_branch, Node* a_node, Node** a_newNode)
{
  ASSERT(a_branch);
  ASSERT(a_node);

  if(a_node->m_count < MAXNODES)  // Split won't be necessary
  {
    a_node->m_branch[a_node->m_count] = *a_branch;
    ++a_node->m_count;
    return false;
  }
  else
  {
    ASSERT(a_newNode);
    
    SplitNode(a_node, a_branch, a_newNode);
	
    return true;
  }
}


// Disconnect a dependent node.
// Caller must return (or stop using iteration index) after this as count has changed
RTREE_TEMPLATE
void RTREE_QUAL::DisconnectBranch(Node* a_node, int a_index)
{
  ASSERT(a_node && (a_index >= 0) && (a_index < MAXNODES));
  ASSERT(a_node->m_count > 0);
  
  // Remove element by swapping with the last element to prevent gaps in array
  a_node->m_branch[a_index] = a_node->m_branch[a_node->m_count - 1];
  
  --a_node->m_count;
}


// Pick a branch.  Pick the one that will need the smallest increase
// in area to accomodate the new rectangle.  This will result in the
// least total area for the covering rectangles in the current node.
// In case of a tie, pick the one which was smaller before, to get
// the best resolution when searching.
RTREE_TEMPLATE
int RTREE_QUAL::PickBranch(Rect* a_rect, Node* a_node)
{
  ASSERT(a_rect && a_node);
  
  bool firstTime = true;
  ELEMTYPEREAL increase;
  ELEMTYPEREAL bestIncr = (ELEMTYPEREAL)-1;
  ELEMTYPEREAL area;
  ELEMTYPEREAL bestArea;
  int best;
  Rect tempRect;

  for(int index=0; index < a_node->m_count; ++index)
  {
    Rect* curRect = &a_node->m_branch[index].m_rect;
    area = CalcRectVolume(curRect);
    tempRect = CombineRect(a_rect, curRect);
    increase = CalcRectVolume(&tempRect) - area;
    if((increase < bestIncr) || firstTime)
    {
      best = index;
      bestArea = area;
      bestIncr = increase;
      firstTime = false;
    }
    else if((increase == bestIncr) && (area < bestArea))
    {
      best = index;
      bestArea = area;
      bestIncr = increase;
    }
  }
  return best;
}


// Combine two rectangles into larger one containing both
RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::CombineRect(Rect* a_rectA, Rect* a_rectB)
{
  ASSERT(a_rectA && a_rectB);

  Rect newRect;

  for(int index = 0; index < NUMDIMS; ++index)
  {
    newRect.m_min[index] = Min(a_rectA->m_min[index], a_rectB->m_min[index]);
    newRect.m_max[index] = Max(a_rectA->m_max[index], a_rectB->m_max[index]);
  }

  return newRect;
}



// Split a node.
// Divides the nodes branches and the extra one between two nodes.
// Old node is one of the new ones, and one really new one is created.
// Tries more than one method for choosing a partition, uses best result.
// Make changes to NodeTransTable and Index units list
RTREE_TEMPLATE
void RTREE_QUAL::SplitNode(Node* a_node, Branch* a_branch, Node** a_newNode)
{
  ASSERT(a_node);
  ASSERT(a_branch);

  // Could just use local here, but member or external is faster since it is reused
  PartitionVars localVars;
  PartitionVars* parVars = &localVars;
  int level;
  // Load all the branches into a buffer, initialize old node
  level = a_node->m_level;
  GetBranches(a_node, a_branch, parVars);                             

  // Find partition
  ChoosePartition(parVars, MINNODES);

  // Put branches from buffer into 2 nodes according to chosen partition
  *a_newNode = AllocNode();
  (*a_newNode)->m_level = a_node->m_level = level;
  (*a_newNode)->unique_id = ++UniqueId;
   
  LoadNodes(a_node, *a_newNode, parVars);

  if(a_node->IsLeaf())
  {
	 mfile.open("savefile", ios::in|ios:: out|std::ios::binary); 
	 int IUtotalsize = 5*sizeof(ELEMTYPE)+sizeof(int)+sizeof(char);
	 char buffer[PAGESIZE];
	 int flag = 0;
	 char op_flag;
	 memset(buffer,'0',PAGESIZE);
	 int_2 table ;
	 std::list<int> splitNodes;
	 splitNodes.push_back(a_node->unique_id);
	 splitNodes.push_back((*a_newNode)->unique_id);
     std::list<int> split_it;
	 update_flag = 1;
	 int j_count=0;
	 ELEMTYPE load_temp[4];
	 int load_id;

	pair<multimap <int,int_2>::iterator, multimap <int,int_2>::iterator > map_it ;//= map_it.second;

	for(int index=0; index < parVars->m_total; ++index)
	{
		
	map_it = NodeTransTable.equal_range(a_node->unique_id);
		
	if(parVars->m_partition[index] == 1)
	{
			
		flag=0;		
		for (std::multimap <int, int_2>::iterator it2 = map_it.first; it2 != map_it.second;)
		{
		
			int seek = (*it2).second.i_0*PAGESIZE;
			mfile.seekg(seek,std::ios_base::beg); 
			mfile.read ((char*)buffer, sizeof(buffer));
	
			memcpy(&load_id,buffer+((*it2).second.i_1* IUtotalsize),sizeof(int));
			memcpy(&load_temp[0],buffer+((*it2).second.i_1* IUtotalsize)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&load_temp[1],buffer+((*it2).second.i_1* IUtotalsize)+sizeof(ELEMTYPE)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&load_temp[2],buffer+((*it2).second.i_1* IUtotalsize)+2*sizeof(ELEMTYPE)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&load_temp[3],buffer+((*it2).second.i_1* IUtotalsize)+3*sizeof(ELEMTYPE)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&op_flag,buffer+((*it2).second.i_1* IUtotalsize)+4*sizeof(ELEMTYPE)+sizeof(int),sizeof(char));
		 
			if((load_temp[0] == parVars->m_branchBuf[index].m_rect.m_min[0])&&(load_temp[1] == parVars->m_branchBuf[index].m_rect.m_min[1])
				&&(load_temp[2] == parVars->m_branchBuf[index].m_rect.m_max[0])&&(load_temp[3] == parVars->m_branchBuf[index].m_rect.m_max[1])&&(load_id==parVars->m_branchBuf[index].m_data))
			{
			
				table.i_0 = (*it2).second.i_0;
				table.i_1 = (*it2).second.i_1; 
				
				split_it.remove(index);
			
				std::multimap <int, int_2>::iterator itcopy = it2;
				itcopy++;
					
				NodeTransTable.erase(it2);
				it2=itcopy;
				
				NodeTransTable.insert(std::pair<int,int_2>((*a_newNode)->unique_id,table));
				
				flag=1;
				break;                                                       
				
			}
				
			else ++it2;
		}
	
			
		if(flag==0)
		{ 
			
				int list_count = 1;
				std::list<int_2>::iterator s_it;

				for(std::list<IndexUnit>::iterator list_it=mylist.begin(); list_it!=mylist.end(); ++list_it)
				{
					
					if(list_count == parVars->m_branchBuf[index].iu_count) 
					{
			
						list_it->iu_id = (int)(*a_newNode)->unique_id;
						flag=2;
						break;
					
					}
					list_count++;
				
				}
				
		
					
		}
		if(flag == 0)
		{
		
			nodebuffer[1] = (*a_newNode);
		}
	}
	

		j_count = 0;

  } 

  

   if(mfile)
   {
		mfile.close();
    }
 
  }
  
  ASSERT((a_node->m_count + (*a_newNode)->m_count) == parVars->m_total);
}

// Calculate the n-dimensional volume of a rectangle
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_QUAL::RectVolume(Rect* a_rect)
{
  ASSERT(a_rect);
  
  ELEMTYPEREAL volume = (ELEMTYPEREAL)1;

  for(int index=0; index<NUMDIMS; ++index)
  {
    volume *= a_rect->m_max[index] - a_rect->m_min[index];
  }
  
  ASSERT(volume >= (ELEMTYPEREAL)0);
  
  return volume;
}


// The exact volume of the bounding sphere for the given Rect
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_QUAL::RectSphericalVolume(Rect* a_rect)
{
  ASSERT(a_rect);
   
  ELEMTYPEREAL sumOfSquares = (ELEMTYPEREAL)0;
  ELEMTYPEREAL radius;

  for(int index=0; index < NUMDIMS; ++index) 
  {
    ELEMTYPEREAL halfExtent = ((ELEMTYPEREAL)a_rect->m_max[index] - (ELEMTYPEREAL)a_rect->m_min[index]) * 0.5f;
    sumOfSquares += halfExtent * halfExtent;
  }

  radius = (ELEMTYPEREAL)sqrt(sumOfSquares);
  
  // Pow maybe slow, so test for common dims like 2,3 and just use x*x, x*x*x.
  if(NUMDIMS == 3)
  {
    return (radius * radius * radius * m_unitSphereVolume);
  }
  else if(NUMDIMS == 2)
  {
    return (radius * radius * m_unitSphereVolume);
  }
  else
  {
    return (ELEMTYPEREAL)(pow(radius, NUMDIMS) * m_unitSphereVolume);
  }
}


// Use one of the methods to calculate retangle volume
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_QUAL::CalcRectVolume(Rect* a_rect)
{
#ifdef RTREE_USE_SPHERICAL_VOLUME
  return RectSphericalVolume(a_rect); // Slower but helps certain merge cases
#else // RTREE_USE_SPHERICAL_VOLUME
  return RectVolume(a_rect); // Faster but can cause poor merges
#endif // RTREE_USE_SPHERICAL_VOLUME  
}


// Load branch buffer with branches from full node plus the extra branch.
RTREE_TEMPLATE
void RTREE_QUAL::GetBranches(Node* a_node, Branch* a_branch, PartitionVars* a_parVars)
{
  ASSERT(a_node);
  ASSERT(a_branch);

  ASSERT(a_node->m_count == MAXNODES);
  
  if(a_node->IsLeaf())
  {
	Node* tmp_node = AllocNode();
	ConstractNode(&tmp_node, a_node->unique_id);
	int index = 0;
	
	while(index < (int)(tmp_node)->m_count)
	{
		a_parVars->m_branchBuf[index] = tmp_node->m_branch[index];
		index++;
	}

/*
	  for(int index=0; index < MAXNODES; ++index)
	{    rtreeloads++;
		a_parVars->m_branchBuf[index] = a_node->m_branch[index];
		printf("data = %d\n",a_parVars->m_branchBuf[index].m_data);
	}*/
  }
  
  else
  {
  
  // Load the branch buffer 
	for(int index=0; index < MAXNODES; ++index)
	{
		a_parVars->m_branchBuf[index] = a_node->m_branch[index];
	}
  }
 
  a_parVars->m_branchBuf[MAXNODES] = *a_branch;
  a_parVars->m_branchCount = MAXNODES + 1;

  // Calculate rect containing all in the set
  a_parVars->m_coverSplit = a_parVars->m_branchBuf[0].m_rect;
  for(int index=1; index < MAXNODES+1; ++index)
  {
    a_parVars->m_coverSplit = CombineRect(&a_parVars->m_coverSplit, &a_parVars->m_branchBuf[index].m_rect);
  }
  a_parVars->m_coverSplitArea = CalcRectVolume(&a_parVars->m_coverSplit);

  InitNode(a_node);
}


// Method #0 for choosing a partition:
// As the seeds for the two groups, pick the two rects that would waste the
// most area if covered by a single rectangle, i.e. evidently the worst pair
// to have in the same group.
// Of the remaining, one at a time is chosen to be put in one of the two groups.
// The one chosen is the one with the greatest difference in area expansion
// depending on which group - the rect most strongly attracted to one group
// and repelled from the other.
// If one group gets too full (more would force other group to violate min
// fill requirement) then other group gets the rest.
// These last are the ones that can go in either group most easily.
RTREE_TEMPLATE
void RTREE_QUAL::ChoosePartition(PartitionVars* a_parVars, int a_minFill)
{
  ASSERT(a_parVars);
  
  ELEMTYPEREAL biggestDiff;
  int group, chosen, betterGroup;
  
  InitParVars(a_parVars, a_parVars->m_branchCount, a_minFill);
  PickSeeds(a_parVars);

  while (((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total)
       && (a_parVars->m_count[0] < (a_parVars->m_total - a_parVars->m_minFill))
       && (a_parVars->m_count[1] < (a_parVars->m_total - a_parVars->m_minFill)))
  {
    biggestDiff = (ELEMTYPEREAL) -1;
    for(int index=0; index<a_parVars->m_total; ++index)
    {
      if(!a_parVars->m_taken[index])
      {
        Rect* curRect = &a_parVars->m_branchBuf[index].m_rect;
        Rect rect0 = CombineRect(curRect, &a_parVars->m_cover[0]);
        Rect rect1 = CombineRect(curRect, &a_parVars->m_cover[1]);
        ELEMTYPEREAL growth0 = CalcRectVolume(&rect0) - a_parVars->m_area[0];
        ELEMTYPEREAL growth1 = CalcRectVolume(&rect1) - a_parVars->m_area[1];
        ELEMTYPEREAL diff = growth1 - growth0;
        if(diff >= 0)
        {
          group = 0;
        }
        else
        {
          group = 1;
          diff = -diff;
        }

        if(diff > biggestDiff)
        {
          biggestDiff = diff;
          chosen = index;
          betterGroup = group;
        }
        else if((diff == biggestDiff) && (a_parVars->m_count[group] < a_parVars->m_count[betterGroup]))
        {
          chosen = index;
          betterGroup = group;
        }
      }
    }
    Classify(chosen, betterGroup, a_parVars);
  }

  // If one group too full, put remaining rects in the other
  if((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total)
  {
    if(a_parVars->m_count[0] >= a_parVars->m_total - a_parVars->m_minFill)
    {
      group = 1;
    }
    else
    {
      group = 0;
    }
    for(int index=0; index<a_parVars->m_total; ++index)
    {
      if(!a_parVars->m_taken[index])
      {
        Classify(index, group, a_parVars);
      }
    }
  }

  ASSERT((a_parVars->m_count[0] + a_parVars->m_count[1]) == a_parVars->m_total);
  ASSERT((a_parVars->m_count[0] >= a_parVars->m_minFill) && 
        (a_parVars->m_count[1] >= a_parVars->m_minFill));
}


// Copy branches from the buffer into two nodes according to the partition.
RTREE_TEMPLATE
void RTREE_QUAL::LoadNodes(Node* a_nodeA, Node* a_nodeB, PartitionVars* a_parVars)
{
  ASSERT(a_nodeA);
  ASSERT(a_nodeB);
  ASSERT(a_parVars);

  for(int index=0; index < a_parVars->m_total; ++index)
   {
    ASSERT(a_parVars->m_partition[index] == 0 || a_parVars->m_partition[index] == 1);
    
    if(a_parVars->m_partition[index] == 0)
    {
      AddBranch(&a_parVars->m_branchBuf[index], a_nodeA, NULL);
    }
    else if(a_parVars->m_partition[index] == 1)
    {
		
      AddBranch(&a_parVars->m_branchBuf[index], a_nodeB, NULL);
    }
  }

}


// Initialize a PartitionVars structure.
RTREE_TEMPLATE
void RTREE_QUAL::InitParVars(PartitionVars* a_parVars, int a_maxRects, int a_minFill)
{
  ASSERT(a_parVars);

  a_parVars->m_count[0] = a_parVars->m_count[1] = 0;
  a_parVars->m_area[0] = a_parVars->m_area[1] = (ELEMTYPEREAL)0;
  a_parVars->m_total = a_maxRects;
  a_parVars->m_minFill = a_minFill;
  for(int index=0; index < a_maxRects; ++index)
  {
    a_parVars->m_taken[index] = false;
    a_parVars->m_partition[index] = -1;
  }
}


RTREE_TEMPLATE
void RTREE_QUAL::PickSeeds(PartitionVars* a_parVars)
{
  int seed0, seed1;
  ELEMTYPEREAL worst, waste;
  ELEMTYPEREAL area[MAXNODES+1];

  for(int index=0; index<a_parVars->m_total; ++index)
  {
    area[index] = CalcRectVolume(&a_parVars->m_branchBuf[index].m_rect);
  }

  worst = -a_parVars->m_coverSplitArea - 1;
  for(int indexA=0; indexA < a_parVars->m_total-1; ++indexA)
  {
    for(int indexB = indexA+1; indexB < a_parVars->m_total; ++indexB)
    {
      Rect oneRect = CombineRect(&a_parVars->m_branchBuf[indexA].m_rect, &a_parVars->m_branchBuf[indexB].m_rect);
      waste = CalcRectVolume(&oneRect) - area[indexA] - area[indexB];
      if(waste > worst)
      {
        worst = waste;
        seed0 = indexA;
        seed1 = indexB;
      }
    }
  }
  Classify(seed0, 0, a_parVars);
  Classify(seed1, 1, a_parVars);
}


// Put a branch in one of the groups.
RTREE_TEMPLATE
void RTREE_QUAL::Classify(int a_index, int a_group, PartitionVars* a_parVars)
{
  ASSERT(a_parVars);
  ASSERT(!a_parVars->m_taken[a_index]);

  a_parVars->m_partition[a_index] = a_group;
  a_parVars->m_taken[a_index] = true;

  if (a_parVars->m_count[a_group] == 0)
  {
    a_parVars->m_cover[a_group] = a_parVars->m_branchBuf[a_index].m_rect;
  }
  else
  {
    a_parVars->m_cover[a_group] = CombineRect(&a_parVars->m_branchBuf[a_index].m_rect, &a_parVars->m_cover[a_group]);
  }
  a_parVars->m_area[a_group] = CalcRectVolume(&a_parVars->m_cover[a_group]);
  ++a_parVars->m_count[a_group]; 
}


// Delete a data rectangle from an index structure.
// Pass in a pointer to a Rect, the tid of the record, ptr to ptr to root node.
// Returns 1 if record not found, 0 if success.
// RemoveRect provides for eliminating the root.
RTREE_TEMPLATE
bool RTREE_QUAL::RemoveRect(Rect* a_rect, const DATATYPE& a_id, Node** a_root)
{
  ASSERT(a_rect && a_root);
  ASSERT(*a_root);

  Node* tempNode;
  ListNode* reInsertList = NULL;

  if(!RemoveRectRec(a_rect, a_id, *a_root, &reInsertList))
  {
    // Found and deleted a data item
    // Reinsert any branches from eliminated nodes
    while(reInsertList)
    {
      tempNode = reInsertList->m_node;
	  
      for(int index = 0; index < tempNode->m_count; ++index)
      {
        InsertRect(&(tempNode->m_branch[index].m_rect),
                   tempNode->m_branch[index].m_data,
                   a_root,
                   tempNode->m_level,tempNode->m_branch[index].iu_count, tempNode->unique_id);
		
      }
      
      ListNode* remLNode = reInsertList;
      reInsertList = reInsertList->m_next;
      FreeNode(remLNode->m_node);
      FreeListNode(remLNode);
    }

    // Check for redundant root (not leaf, 1 child) and eliminate
    if((*a_root)->m_count == 1 && (*a_root)->IsInternalNode())
    {	
      tempNode = (*a_root)->m_branch[0].m_child;
      
      ASSERT(tempNode);
      FreeNode(*a_root);
      *a_root = tempNode;
    }
    return false;
  }
  else
  {
    return true;
  }
}


// Delete a rectangle from non-root part of an index structure.
// Called by RemoveRect.  Descends tree recursively,
// merges branches on the way back up.
// Returns 1 if record not found, 0 if success.
RTREE_TEMPLATE
bool RTREE_QUAL::RemoveRectRec(Rect* a_rect, const DATATYPE& a_id, Node* a_node, ListNode** a_listNode)
{
	

  ASSERT(a_rect && a_node && a_listNode);
  ASSERT(a_node->m_level >= 0);

  if(a_node->IsInternalNode())  // not a leaf node
  {
	  
    for(int index = 0; index < a_node->m_count; ++index)
    {
      if(Overlap(a_rect, &(a_node->m_branch[index].m_rect)))
      {
        if(!RemoveRectRec(a_rect, a_id, a_node->m_branch[index].m_child, a_listNode))
        {
          if(a_node->m_branch[index].m_child->m_count >= MINNODES)
          {
            // child removed, just resize parent rect
            a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
          }
          else
          {
            // child removed, not enough entries in node, eliminate node
            ReInsert(a_node->m_branch[index].m_child, a_listNode);
	
            DisconnectBranch(a_node, index); // Must return after this call as count has changed
          }
          return false;
        }
      }
    }
    return true;
  }
  else // A leaf node
  {

	  if(a_node->unique_id < 0)
	  {
	  a_node->unique_id = ++UniqueId;    
	  
	  }
	    
	 nodebuffer[1] = a_node;
	
    for(int index = 0; index < a_node->m_count; ++index)
    {

      if(a_node->m_branch[index].m_child == (Node*)a_id)
      {
		  IUcount++;
		   nodebuffer[1] = a_node;
        DisconnectBranch(a_node, index); // Must return after this call as count has changed
		 
        return false;
      }
	  
    }

    return true;
  }
}


// Decide whether two rectangles overlap.
RTREE_TEMPLATE
bool RTREE_QUAL::Overlap(Rect* a_rectA, Rect* a_rectB)
{
  ASSERT(a_rectA && a_rectB);

  for(int index=0; index < NUMDIMS; ++index)
  {
    if (a_rectA->m_min[index] > a_rectB->m_max[index] ||
        a_rectB->m_min[index] > a_rectA->m_max[index])
    {
      return false;
    }
  }
  return true;
}


// Add a node to the reinsertion list.  All its branches will later
// be reinserted into the index structure.
RTREE_TEMPLATE
void RTREE_QUAL::ReInsert(Node* a_node, ListNode** a_listNode)
{
  ListNode* newListNode;
 // printf("reinsert %d\n",a_node->m_branch->iu_count);
  newListNode = AllocListNode();
  newListNode->m_node = a_node;
  newListNode->m_next = *a_listNode;
  *a_listNode = newListNode;
}


// Search in an index tree or subtree for all data retangles that overlap the argument rectangle.
RTREE_TEMPLATE
bool RTREE_QUAL::Search(Node* a_node, Rect* a_rect, int& a_foundCount, bool __cdecl a_resultCallback(DATATYPE a_data, void* a_context), void* a_context)
{
  ASSERT(a_node);
  ASSERT(a_node->m_level >= 0);
  ASSERT(a_rect);

  if(a_node->IsInternalNode()) // This is an internal node in the tree
  {
    for(int index=0; index < a_node->m_count; ++index)
    {
      if(Overlap(a_rect, &a_node->m_branch[index].m_rect))
      {
        if(!Search(a_node->m_branch[index].m_child, a_rect, a_foundCount, a_resultCallback, a_context))
        {
          return false; // Don't continue searching    
        }
      }
    }
  }
  else // This is a leaf node
  {

	Node* tmp_node;
	ConstractNode(&tmp_node, a_node->unique_id);

    for(int index=0; index < tmp_node->m_count; ++index)
    {
      if(Overlap(a_rect, &tmp_node->m_branch[index].m_rect))
      {
        DATATYPE& id = tmp_node->m_branch[index].m_data;
        
        // NOTE: There are different ways to return results.  Here's where to modify
        if(&a_resultCallback)
        {
          ++a_foundCount;
          if(!a_resultCallback(id, a_context))
          {
            return false; // Don't continue searching
          }
        }
      }
    }
  }

  return true; // Continue searching
}
RTREE_TEMPLATE
void RTREE_QUAL::ConstractNode(Node **a_node, int node_id)
{
	 mfile.open("savefile",ios::out|ios::in|std::ios::binary);  
	   
	 if(!mfile)
	   {
	   mfile.clear();
	   }

	   *a_node = AllocNode();
	   (*a_node)->m_level = 0;
	   (*a_node)->m_count = 0;
	   (*a_node)->unique_id = node_id;
	   char buffer[PAGESIZE];

	   struct tempbranch
	   {           
		   int t_id;
		   ELEMTYPE t_min[2];
		   ELEMTYPE t_max[2];
		   char t_flag;
		   int t_iucount;
	   };

	   std::list< tempbranch> mylist_i;   //list for inserts
	   std::list< tempbranch> mylist_r;   //list for removes

	   tempbranch ins;
	   tempbranch rem;
	   int i=1;
	   int IUtotalsize = 5*sizeof(ELEMTYPE)+sizeof(int)+sizeof(char);
	   int flag = 0;
	   memset(buffer,'0',PAGESIZE);
  
	   pair<multimap <int,int_2>::iterator, multimap <int,int_2>::iterator > map_it ;
	   map_it = NodeTransTable.equal_range(node_id);

      int index = 0;
      char op;
	 
      for (multimap <int, int_2>::iterator it2 = map_it.first; it2 != map_it.second; ++it2)
	  {
		int seek = (*it2).second.i_0*PAGESIZE;

		mfile.seekg(seek,std::ios_base::beg); 
		mfile.read ((char*)buffer, sizeof(buffer));
		
		if (mfile.fail())
		{
            mfile.clear();
			mfile.seekg(seek,std::ios_base::beg);
            mfile.read ((char*)buffer, sizeof(buffer));
         }

		flag=1;
		int newid=-1;
	
		memcpy(&op,buffer+((*it2).second.i_1* IUtotalsize)+5*sizeof(int),sizeof(char));
	
		if(op=='i')
		{
			memcpy(&ins.t_id,buffer+((*it2).second.i_1* IUtotalsize),sizeof(int));
			memcpy(&ins.t_min[0],buffer+((*it2).second.i_1* IUtotalsize)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&ins.t_min[1],buffer+((*it2).second.i_1* IUtotalsize)+sizeof(ELEMTYPE)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&ins.t_max[0],buffer+((*it2).second.i_1* IUtotalsize)+2*sizeof(ELEMTYPE)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&ins.t_max[1],buffer+((*it2).second.i_1* IUtotalsize)+3*sizeof(ELEMTYPE)+sizeof(int),sizeof(ELEMTYPE));
			ins.t_flag = op;
			mylist_i.push_back(ins);
		}

		if(op=='r')
		{
			memcpy(&rem.t_id,buffer+((*it2).second.i_1* IUtotalsize),sizeof(int));
			memcpy(&rem.t_min[0],buffer+((*it2).second.i_1* IUtotalsize)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&rem.t_min[1],buffer+((*it2).second.i_1* IUtotalsize)+sizeof(ELEMTYPE)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&rem.t_max[0],buffer+((*it2).second.i_1* IUtotalsize)+2*sizeof(ELEMTYPE)+sizeof(int),sizeof(ELEMTYPE));
			memcpy(&rem.t_max[1],buffer+((*it2).second.i_1* IUtotalsize)+3*sizeof(ELEMTYPE)+sizeof(int),sizeof(ELEMTYPE));
			rem.t_flag = op;
			mylist_r.push_back(rem);
		
		}

		index++;
	
	  }
	
	  for(std::list<IndexUnit>::iterator list_it=mylist.begin(); list_it!=mylist.end(); ++list_it)
	  {
	  
		  if(list_it->IUNode->unique_id ==  (*a_node)->unique_id)
		  {
			  if((char)list_it->op_flag == 'i')
			  {
				  ins.t_id = (int)list_it->id;
				  ins.t_min[0] = list_it->minimal_bounding_box.m_min[0];
				  ins.t_min[1] = list_it->minimal_bounding_box.m_min[1];
				  ins.t_max[0] = list_it->minimal_bounding_box.m_max[0];
				  ins.t_max[1] = list_it->minimal_bounding_box.m_max[1];
				  ins.t_flag = 'i';
				  ins.t_iucount = i;
				  mylist_i.push_back(ins);
				
			  }

			  else if((char)list_it->op_flag == 'r')
			  {
				  rem.t_id = (int)list_it->id;
				  rem.t_min[0] = list_it->minimal_bounding_box.m_min[0];
				  rem.t_min[1] = list_it->minimal_bounding_box.m_min[1];
				  rem.t_max[0] = list_it->minimal_bounding_box.m_max[0];
				  rem.t_max[1] = list_it->minimal_bounding_box.m_max[1];
				  rem.t_flag = 'r';
				  rem.t_iucount = i;
				  mylist_r.push_back(rem);
			  
			  }
		  }
	  
	  ++i;
	  }

	 std::list<tempbranch>::iterator i_it;
	 std::list<tempbranch>::iterator r_it;

	for(std::list<tempbranch>::iterator r_it=mylist_r.begin(); r_it!=mylist_r.end(); r_it++){

	
		for(std::list<tempbranch>::iterator i_it=mylist_i.begin(); i_it!=mylist_i.end(); ++i_it)
		{
			if((r_it->t_id == i_it->t_id)&&(r_it->t_min[0]==i_it->t_min[0])&&(r_it->t_min[1]==i_it->t_min[1])
				&&(r_it->t_max[0]==i_it->t_max[0])&&(r_it->t_max[1]==i_it->t_max[1]))
			{ 
				std::list<tempbranch>::iterator itcopy = i_it;
			
				itcopy++;
				mylist_i.erase(i_it);
				i_it=itcopy;
				break;			
			}
		}
	 }

	index = 0;
	
	for(std::list<tempbranch>::iterator i_it=mylist_i.begin(); i_it!=mylist_i.end(); ++i_it)
	{

		(*a_node)->m_branch[index].iu_count = i_it->t_iucount;
		(*a_node)->m_branch[index].m_data = i_it->t_id;
		(*a_node)->m_branch[index].m_rect.m_min[0] = i_it->t_min[0];
		(*a_node)->m_branch[index].m_rect.m_min[1] = i_it->t_min[1];
		(*a_node)->m_branch[index].m_rect.m_max[0] = i_it->t_max[0];
		(*a_node)->m_branch[index].m_rect.m_max[1] = i_it->t_max[1];
		(*a_node)->m_count++;
		index++;
	}
	if(mfile)
    {
		mfile.close();
    }
	
}

#undef RTREE_TEMPLATE
#undef RTREE_QUAL

#endif //RTREE_H
