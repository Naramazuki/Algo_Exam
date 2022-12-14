""" TP3 where we code a global alignment program
    based on Needleman-Wunsch algorithm """





def simple_display(seq_top, seq_left, score):
    """ Do a simple display """
    print("Seq_top:  {}\nSeq_left: {}\nScore: {}".format(seq_top[::-1], seq_left[::-1], score))

def nice_display(seq_top, seq_left, score):
    """ Do a nice display """
    # What will be printed
    to_print = ""

    # Print seq_top
    # For each positions in a sequence
    for i, _ in enumerate(seq_top):
        # Add the corresponding letter (reverse order)
        to_print += seq_top[-(i+1)]
    # Print a next line after printing the seq_top
    to_print += "\n"

    # Print middle line
    # For each positions in a sequence
    for i, _ in enumerate(seq_top):
        # If it is a match between the two sequences
        if seq_top[-(i+1)] == seq_left[-(i+1)]:
            # Print a pipe
            to_print += "|"
        # Otherwise
        else:
            # Print a space
            to_print += " "
    # Print a next line after printing the middle line
    to_print += "\n"

    # Print seq_left
    # For each positions in a sequence
    for i, _ in enumerate(seq_top):
        # Add the corresponding letter (reverse order)
        to_print += seq_left[-(i+1)]

    # Add the score at the end
    to_print += "\nScore: {}\n".format(score)
    # Print everything!
    print(to_print)

class Cell():
    def __init__(self, score=None, prev_pos=None):
        self.score = score
        self.prev_pos = prev_pos




class DynamicMatrix:
    """ Class to generate an empty matrix """
    def __init__(self, seq_top, seq_left, match, mismatch, indel):
        # Init all "self" variables
        self.seq_top = seq_top
        self.seq_left = seq_left
        self.match = match
        self.mismatch = mismatch
        self.indel = indel
        # Create the matrix of Cell()
        self.matrix=[]        
        for W in range(len(seq_left)+1):
            self.matrix.append([])
            liste=[]
            for j in range(len(seq_top)+1):
                c=Cell()
                liste.append(c)
            
            self.matrix[W]=liste
        

    # self representation for print
    def __repr__(self):
        # What will be returned
        return "Scores:\n{}\nPrev_pos:\n{}\n\n".format(self.print_scores(), self.print_prev_pos())
    #compare two nucleotides
    def compare(self,ntd_A,ntd_B):
        if(ntd_A == ntd_B):
            return self.match
        return self.mismatch
    # self representation for print
    def print_scores(self):
        """ Output the values of the matrix """
        # What will be returned
        ret_scores = ".  . "
        # Print top_seq
        for i in self.seq_top:
            ret_scores += "  {} ".format(i)
        # New line
        ret_scores += "\n"
        # For each line
        for ind, i in enumerate(self.matrix):
            # Print seq_left
            if ind > 0:
                ret_scores += "{} ".format(self.seq_left[ind-1])
            else:
                ret_scores += ". "
            # For each column
            for j in i:
                # If this cell has no value
                if j.score is None:
                    # Add a dot to the return
                    ret_scores += (" . ")
                # If this cell is not empty
                else:
                    # Add its content to the return
                    tmp_val = str(j.score)
                    if len(tmp_val) == 1:
                        ret_scores += " " + tmp_val + " "
                    if len(tmp_val) == 2:
                        ret_scores += tmp_val + " "
                    if len(tmp_val) == 3:
                        ret_scores += tmp_val
                # Always add a space after the value we add
                ret_scores += " "
            # End of this line, go to next line
            ret_scores += "\n"
        # Return the content of the Matrix
        return ret_scores

    # self representation for print
    def print_prev_pos(self):
        """ Output the values of the matrix """
        # What will be returned
        ret_prev_pos = ".   .  "
        # Print top_seq
        for i in self.seq_top:
            ret_prev_pos += "     {} ".format(i)
        # New line
        ret_prev_pos += "\n"
        # For each line
        for ind, i in enumerate(self.matrix):
            # Print seq_left
            if ind > 0:
                ret_prev_pos += "{} ".format(self.seq_left[ind-1])
            else:
                ret_prev_pos += ".   "
            # For each column
            for j in i:
                # If this cell has no value
                if j.prev_pos is None:
                    # Add a dot to the return
                    ret_prev_pos += (".   ")
                # If this cell is not empty
                else:
                    # Add its content to the return
                    tmp_val = str(j.prev_pos)
                    ret_prev_pos += tmp_val
                # Always add a space after the value we add
                ret_prev_pos += " "
            # End of this line, go to next line
            ret_prev_pos += "\n"
        # Return the content of the Matrix
        return ret_prev_pos

    def initialize(self):
        """ Initialize the matrix, i.e. fill the first line and column """
        # First cell is 0
        self.matrix[0][0].score=0
        self.matrix[0][0].prev_pos=[]


        for i in range(1,len(self.seq_top)+1):
            self.matrix[0][i].score=self.matrix[0][i-1].score+self.indel
            self.matrix[0][i].prev_pos=[0,i-1]
        for j in range(1,len(self.seq_left)+1):
            self.matrix[j][0].score=self.matrix[j-1][0].score+self.indel
            self.matrix[j][0].prev_pos=[j,0]

        
    def fill_matrix(self):
        """ Fill-up the matrix """
        for i in range(1,len(self.seq_left)+1):
            for j in range(1,len(self.seq_top)+1):
                self.matrix[i][j].prev_pos=[]

                diag=self.matrix[i-1][j-1].score+self.compare(self.seq_top[j-1],self.seq_left[i-1])
                gap_left=self.matrix[i][j-1].score+self.indel
                gap_top=self.matrix[i-1][j].score+self.indel
                maxi=max( diag,gap_top,gap_left)
                self.matrix[i][j].score=maxi
                if maxi==gap_top:
                    self.matrix[i][j].prev_pos.append([i-1,j])
                if maxi==gap_left:
                    self.matrix[i][j].prev_pos.append([i,j-1])
                if maxi==diag:
                    self.matrix[i][j].prev_pos.append([i-1,j-1])
    def global_alignment2(self,j,i,listetop,listeleft):
            """ Make a global alignment of two sequences """
            
            intstop=j*i
            
            top=""
            left=""
            d=len(self.matrix[j][i].prev_pos)
        

            
        
        
            

            if i>0 and j>0 and intstop>0:
                for n in range (d):
                    v=self.matrix[j][i].prev_pos[n][0]
                    w=self.matrix[j][i].prev_pos[n][1]
                    print(n)
                    if(w==i-1 and v==j-1):
                        listeleft+=self.seq_left[j-1]
                        listetop+=self.seq_top[i-1]
                        i-=1
                        j-=1
                        top,left=self.global_alignment2(j,i,listetop,listeleft)
                        
                        listetop=listetop+top
                        listeleft=listeleft+left
            
                    elif(w==i-1 and v==j):
                        listeleft+="-"
                        listetop+=self.seq_top[i-1]
                        i-=1
                        top,left=self.global_alignment2(j,i,listetop,listeleft)
                        
                        listetop=listetop+top
                        listeleft=listeleft+left

                        
                    elif(w==i and v==j-1):
                        listeleft+=self.seq_left[j-1]
                        listetop+="-"
                        j-=1
                        top,left=self.global_alignment2(j,i,listetop,listeleft)
                        
                        listetop=listetop+top
                        listeleft=listeleft+left

            intstop-=1
       # faire partie en iif qui va faire appel a deux focntion différentes selon si on se trouve dans le cas d'un embranchement ou danns le cas normal.
       # argurment de la fonction sera la position en coordonnées de l'embranchement 

            return(listetop,listeleft)

    def global_alignment(self):
        """ Make a global alignment of two sequences """
        j,i=len(self.seq_left),len(self.seq_top)
        
        
        seq_top=[]
        seq_left=[]
        score=self.matrix[j][i].score
        listetop=""
        listeleft=""
        seq_top,seq_left=self.global_alignment2(j,i,listetop,listeleft)
       
        
        return(seq_top,seq_left,score)



            


def main():
    """ The main of TP3"""
    #mat = DynamicMatrix("ACGGCTATTCA", "ACTGTAGGGT", 2, -1, -2)
    mat = DynamicMatrix("ACGGCTAT", "ACTGTAG", 2, -1, -2)

    mat.initialize()
    mat.fill_matrix()
    
    print(mat)
   


    al_seq_top, al_seq_left, score = mat.global_alignment()
    #nice_display(al_seq_top, al_seq_left, score)

# Launch the main
main()
# Exit without error
exit(0)
# Always put one extra return line
