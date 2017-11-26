#include "../CollocationPoint.hpp"

//------- create data structure for a node ---------//----------------------------------//
struct node {                                       //                                  //
    CollocationPoint collPnt;                       // the collocation point as object  //
    node* left;                                     // left child node                  //
    node* right;                                    // right child node                 //

};

//------- function to create a new node ------------//----------------------------------//
node* get_new_node(CollocationPoint collPnt){       //                                  //
    node* new_node = new node();                    // create pointer to new node       //
    new_node->collPnt = collPnt;                    // set the data for the new node    //
    new_node->left = NULL;                          // initialize left child as empty   //
    new_node->right = NULL;                         // initialize right child as empty  //
    return new_node;                                // return the new node              //

};

//------- insert data into the next node -----------//----------------------------------//
node* insert(node* root, double data) {             //                                  //
    if (root==NULL) {                               // if root is empty, give structure //
        root = get_new_node(data);                  //                                  //
    } else if (data <= root->data) {                // 
        root->left = insert(root->left,data);
    } else if (data > root->data) {
        root->right = insert(root->right,data);
    }
    return root;
};

void print_tree(struct node* node) {
    if (node==NULL) return;
    print_tree(node->left);
    print_tree(node->right);
    cout << node->data << endl;
 //    printf("%f",node->data);
};

bool lookup(node* root,int data) {
    if (root == NULL) return false;
    else if (root->data == data) return true;
    else if (data <= root->data) return lookup(root->left,data);
    else return lookup(root->right,data);
};

int main() {
    node* root = NULL;          // setting tree as empty
    root = insert(root, -1);
    root = insert(root,-0.25);
    print_tree(root);
    return 0;
}
    
