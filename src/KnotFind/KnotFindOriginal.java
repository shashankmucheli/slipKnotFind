/*********************************************************************************************************************
 * FAQ
 * Q: What program is needed to run the code?
 * A: The project was built in NetBeans IDE. Download latest version of NetBeans form https://netbeans.org/downloads/
 *    Also, I believe you can import this project by Eclipse. 
 *    Download Eclipse form https://www.eclipse.org/downloads/packages/eclipse-standard-432/keplersr2
 *    
 * Q: Where should I place the PDB file?
 * A: The PDB file should be in the root director of the project. 
 *    One PDB file "1ALK_A.pdb" is provided.
 * 
 * Q: How can I run my own PDB file?
 * A: You can run your own PDB file through modefying the static final String PDB variable inside class SlipknotFind.
 *    Of course, you can also modify the PATH variable to your PDB folder.
 * 
 * Q: Do I need to pre-process PDB file? What do you mean by "1ALK_A"?
 * A: PDB "1ALK" has two chains (chain A and chain B). In this case, you have to seperate the two chains into two files.
 *    See the differences between "1ALK.pdb" and "1ALK_A.pdb" inside the project.
 *    For more information about PDB file. http://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)
 * 
 * For more questions: yzhang5@umassd.edu 
 ***********************************************************************************************************************/


package KnotFind;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import slipknotfind.SlipknotFind;


public class KnotFindOriginal {
   
   
   static final double TOLERANCE=0.0003;
   static List<Triangle> tri=new ArrayList();
   static List<Res> res=new ArrayList();
   static boolean _byArea=false;
   static final String PATH="";     //the path of your pdb file.     e.g. PATH="PDB/";
   static final String PDB="1JS1.pdb";
   int first_atom,last_atom = 0;
   
   
   public static void main(String[] args) {
      KnotFindOriginal sf=new KnotFindOriginal(); 
      sf.initResidual(PDB);
      sf.slipknotFind(res);
   }
   
   void initResidual(String name){
      String fileName=PATH+name;
      File file=new File(fileName);
      BufferedReader reader=null;
      try{
         reader=new BufferedReader(new FileReader(file));
         String temp;
         while((temp=reader.readLine())!= null){
            String[] s=temp.split("\\s+");
            if(s[0].equals("ATOM"))
               if(s[2].equals("CA")){
                   if(s[4].equals("X")){
                       if(first_atom == 0) { first_atom = Integer.parseInt(s[5]); }
                  Res r=new Res();
                  r.index=Integer.parseInt(s[5]);
                  last_atom = Integer.parseInt(s[5]);
                  r.x=Double.parseDouble(s[6]);
                  r.y=Double.parseDouble(s[7]);
                  r.z=Double.parseDouble(s[8]);
                  res.add(r);
                   }
                }
                   
         }
      }catch(IOException e){
         e.printStackTrace();
      }
      System.out.println("residual.size=" + res.size());
//      System.out.println("*****"+res.size()+" CA information**********");
//      for(int i=0;i<res.size();i++){
//         System.out.println(res.get(i).index+":\t"+res.get(i).x+"\t"+res.get(i).y+"\t"+res.get(i).z);
//      }
//      System.out.println("***************************************");
   }
   
   public boolean knotFind(List res){
      simplify(res);
      if(res.size()>2){
         dblCheck(res);
      }
      
      if (res.size() == 2) {
         return false;
      } 
      else {
//         System.out.println("final res.size = " + res.size());
//         System.out.println("the remained chain numbers are:");
//         for (int i = 0; i < res.size(); i++) {
//            System.out.print(res.get(i).index + " ");
//         }
//         System.out.println();
         return true;
      }

   }
   
   public void slipknotFind(List<Res> res){
      List<Res> r;  //a part of res
      boolean _knotInR=false;
      boolean _knotInRes=true;
      int k3=0;
      int k2=0;
      
      //init r
      for(int i=0;i<res.size()-2;i++){
         if(_knotInR) break;
         for(int j=i+2;j<res.size();j++){
            if(_knotInR) break;
            System.out.println("checking residues from " + i + " to "+ j +"  ");
            r=new ArrayList();
            for(int p=0;p<j+1;p++){
               r.add(res.get(p));
            }
            initTriangle(r);

            if(knotFind(r)){
               _knotInR=true;
               k3=i;
               k2=j;
               System.out.println("find a knot between "+ k3 + " and "+k2);
            }
         }
      }
      
      for(int i=k3+1;i<=k2;i++){
         r=new ArrayList();
         for(int p=i;p<=k2;p++){
            r.add(res.get(p));
         }
         initTriangle(r);
         
         System.out.println("checking residues from " + i + " to "+ k2 +"  ");
         if(!knotFind(r)){
            k3=i-1;
            System.out.println("find smallest knot between "+k3 +" to "+ k2);
            break;
         }
      }
      
      
      if(_knotInR){  
         System.out.println("***********************************************");
         System.out.println("now check if the knot could be untied eventually");
         int k1=-1;
         for(int i=k2+1;i<res.size();i++){
            System.out.println("checking residues from "+k3+" to "+ i);
            r=new ArrayList();
            for(int p=k3;p<i;p++){
               r.add(res.get(p));
            }  
            initTriangle(r);
            if(!knotFind(r)){
               _knotInRes=false;
               k1=i;
               break;
            }
         }
         
         if(_knotInRes){      
            for(int i = k3 - 1; i >= 0; i--) {
               System.out.println("checking residues from " + i + " to " + k2);
               r = new ArrayList();
               for (int p = k2; p >= i; p--) {
                  r.add(res.get(p));
               }

               initTriangle(r);
               if (!knotFind(r)) {
                  _knotInRes = false;
                  k1 = i;
                  break;
               }
            }
         }

         if(_knotInR == true && _knotInRes == false){
            System.out.println("find a slipknot: k3="+k3+"  k2="+k2+"  k1="+ k1);
         }
         else{
            System.out.println("this chain only has a knot between "+ k3+ " and "+k2);
         }
      }
      else{
         System.out.println("not knots, no slipknots");
      }
   }
   
   void simplify(List<Res> res){
      int numInTri=0;
      while(res.size()>2){
         if(numInTri>=tri.size()){
            //cannot simplify any more
            ArrayList list = new ArrayList();
            for(Res s : res){
                list.add(s.index);
            }
            System.out.println(list);
            createpdb(list);
            return;
         }

         int res1 = tri.get(numInTri).res1;
         Plane plane = new Plane(res,res1);
         if (!findIntersect(plane,res)) {
            deleteResidual(res, res1 + 1);
            numInTri = 0;
         } else {
            numInTri++;
         }
      }
   }
   
   boolean findIntersect(Plane plane, List<Res> res){
      List<Integer> trouble=new ArrayList();
      
      int point1=0;
      //find all the trouble segments
      while(point1<res.size()-1){
         int point2=point1+1;
         //check if point1 or point2 is part the triangle
         if(point1 == plane.residual || point1 == plane.residual+1 || point1 == plane.residual+2 || point2 == plane.residual){
            point1++; continue;
         }
            
         //check if p1 and p2 are in the same side of the plane
         double p1Distance, p2Distance;
         p1Distance=plane.a*res.get(point1).x + plane.b * res.get(point1).y + plane.c * res.get(point1).z + plane.d;
         p2Distance=plane.a*res.get(point2).x + plane.b * res.get(point2).y + plane.c * res.get(point2).z + plane.d;
         //p1 and p2 are on the same side
         if(p1Distance>TOLERANCE && p2Distance > TOLERANCE || p1Distance < -TOLERANCE && p2Distance < -TOLERANCE){
            point1++; continue;
         }
         //now find a trouble point the starting and ending point are on the different side of the plane
         trouble.add(point1++);
      }
      
      boolean _intersect=false;
      //check if the trouble segments are really intersect the plane
      for(int i=0;i<trouble.size();i++){
         int p1=trouble.get(i);
         int p2=trouble.get(i)+1;
         
         double denom=plane.a * (res.get(p1).x - res.get(p2).x)
                        + plane.b * (res.get(p1).y - res.get(p2).y)
                           + plane.c * (res.get(p1).z - res.get(p2).z);
         
                 
         if(Math.abs(denom)<TOLERANCE){         //denominator is 0.   the segment is inside the plane
            System.out.println("****the segment is inside the plane now****"+ Math.abs(denom));
            //
            if(segmentIntersectTriangle(plane,p1,res)){
               return true;
            }
            if(segmentIntersectTriangle(plane,p2,res)){
               return true;
            }
         }
         else{                                  //the segment crosses the plane
            double mu=( plane.d + plane.a * res.get(p1).x + 
                           plane.b * res.get(p1).y + 
                              plane.c * res.get(p1).z ) / denom;
            
            if( mu < TOLERANCE || mu > 1.0+TOLERANCE){
               //System.out.println("mu > 1 or mu < 0 on plane: "+plane.residual+" point: "+point1);
               continue;
            }
            //cal the cross point
            Res crossPoint=new Res();
            crossPoint.x=res.get(p1).x +(mu * (res.get(p2).x - res.get(p1).x));
            crossPoint.y=res.get(p1).y +(mu * (res.get(p2).y - res.get(p1).y));
            crossPoint.z=res.get(p1).z +(mu * (res.get(p2).z - res.get(p1).z));
            //check if the cross point is inside the triangle
            if(crossPointInsideTriangle(plane,crossPoint,res)){
               //System.out.println("cross point: \t"+crossPoint.x+"\t"+crossPoint.y+"\t"+crossPoint.z);
               return _intersect=true;
            }
         }
      }//end check if the trouble segments are really intersect the plane
      
      return _intersect;
   }
   
   boolean segmentIntersectTriangle(Plane plane, int p1, List<Res> res){
      if(crossPointInsideTriangle(plane,res.get(p1),res))
         return true;
      if(crossPointInsideTriangle(plane,res.get(p1+1),res))
         return true;
      
      return false;
   }
   
   boolean crossPointInsideTriangle(Plane plane, Res p, List<Res> res){
      Res a=res.get(plane.residual);
      Res b=res.get(plane.residual+1);
      Res c=res.get(plane.residual+2);
      
      double ab=Math.pow(b.x-a.x,2) + Math.pow(b.y-a.y,2) + Math.pow(b.z-a.z,2);
      double bc=Math.pow(c.x-b.x,2) + Math.pow(c.y-b.y,2) + Math.pow(c.z-b.z,2);
      double ca=Math.pow(a.x-c.x,2) + Math.pow(a.y-c.y,2) + Math.pow(a.z-c.z,2);
      
      double pa=Math.pow(p.x-a.x,2) + Math.pow(p.y-a.y,2) + Math.pow(p.z-a.z,2);
      double pb=Math.pow(p.x-b.x,2) + Math.pow(p.y-b.y,2) + Math.pow(p.z-b.z,2);
      double pc=Math.pow(p.x-c.x,2) + Math.pow(p.y-c.y,2) + Math.pow(p.z-c.z,2);
      
      double temp;
      
      temp=(pb+pc-bc)/(2*Math.sqrt(pb*pc));
      double bpc=Math.acos(temp);
      
      temp=(pb+pa-ab)/(2*Math.sqrt(pa*pb));
      double apb=Math.acos(temp);
      
      temp=(pa+pc-ca)/(2*Math.sqrt(pa*pc));
      double cpa=Math.acos(temp);
      
      double total=bpc+apb+cpa;
      //System.out.println("total: "+total);
      if(Math.abs(total-2*Math.PI) < TOLERANCE){
         //System.out.println("plane:"+ plane.residual+" has crossings ");
         return true;
      }
      else{
         //System.out.println("outside the triangle***********oh yeah");
         return false;
      }
        
   }
   
   void deleteResidual(List<Res> res, int res2){
      //delete 3 triangles
      int del=0;
      int add=0;
      for(int i=0;i<tri.size();i++){
         if(tri.get(i).res1==res2-2){
            tri.remove(i);del++;
            break;
         }
      }
      for(int i=0;i<tri.size();i++){
         if(tri.get(i).res1==res2-1){
            tri.remove(i);del++;
            break;
         }
      }
      
      for(int i=0;i<tri.size();i++){
         if(tri.get(i).res1==res2){
            tri.remove(i);del++;
            break;
         }
      }
      
      for(int i=0;i<tri.size();i++){
         if(tri.get(i).res1>res2)
            tri.get(i).res1--;
      }
      
      //remove the residual from res chain
      //System.out.println("Deleting...  "+res.get(res2).index);
     // deleted_atoms.add(res.get(res2).index);
      res.remove(res2);
      
      if(res2>1){
      //add 2 new triangles
         Triangle newTri1=new Triangle();
         newTri1.res1=res2-2;
         newTri1.distance=calDistance(res, res2-2);
         insertIntoTri(newTri1); add++;
      }
      if(res2<res.size()-1){
         Triangle newTri2=new Triangle();
         newTri2.res1=res2-1;
         newTri2.distance=calDistance(res, res2-1);
         insertIntoTri(newTri2); add++;
      }
      if(del-add!=1)
         System.out.println("delete "+del+" triangles and add "+ add +" tri");
   } 
   
   void initTriangle(List<Res> res){
      tri.clear();
      int[] res1 = new int[res.size() - 2];
      double[] distance = new double[res.size() - 2];
      
      for(int i=0;i<res.size()-2;i++){
         res1[i]=i;
         distance[i]=calDistance(res, i);
      }
      Triangle t=new Triangle();
      t.res1=res1[0];
      t.distance=distance[0];
      tri.add(t);
      
      for(int i=1;i<res1.length;i++){
         t=new Triangle();
         t.res1=res1[i];
         t.distance=distance[i];
         insertIntoTri(t);
      }
   }
   
  void insertIntoTri(Triangle newTri){
      int index=-1;
      for(int i=tri.size()-1;i>=0;i--){
         if(newTri.distance>tri.get(i).distance){
            index=i;
            break;
         }
      }
      tri.add(index+1,newTri);
   }
   
  double calDistance(List<Res> res, int res1){
      double distance;
      int a=res1;
      int b=a+1;
      int c=a+2;
      if(!_byArea){
         distance=Math.pow((res.get(a).x-res.get(c).x),2)
                  +Math.pow((res.get(a).y-res.get(c).y),2)
                   +Math.pow((res.get(a).z-res.get(c).z),2);
      }
      else{
         
         double ab,bc,ca;
         ab=Math.sqrt(  Math.pow(res.get(b).x-res.get(a).x, 2)
                           +Math.pow(res.get(b).y-res.get(a).y,2)
                              +Math.pow(res.get(b).z-res.get(a).z, 2));
         bc=Math.sqrt(  Math.pow(res.get(b).x-res.get(c).x, 2)
                           +Math.pow(res.get(b).y-res.get(c).y,2)
                              +Math.pow(res.get(b).z-res.get(c).z, 2));
         ca=Math.sqrt(  Math.pow(res.get(c).x-res.get(a).x, 2)
                           +Math.pow(res.get(c).y-res.get(a).y,2)
                              +Math.pow(res.get(c).z-res.get(a).z, 2));
         double s=((ab+bc+ca)/2);
         distance=(s*(s-ab)*(s-bc)*(s-ca));
      }
      
      return distance;
      
   }
   
   void dblCheck(List<Res> res){
//      System.out.println("******************************");
//      System.out.println("now double checking");
      _byArea=true;
      initTriangle(res);
      simplify(res);       
      _byArea=false;
   }

    private void createpdb(ArrayList list) {
        String fileName=PATH+PDB;
      File file=new File(fileName);
      BufferedReader reader=null;
      BufferedReader tmp=null;
      /*System.out.println("k1 value : " + k1);
      System.out.println("k3 value : " + k3);*/
      try{
         String trim_filename = PATH+"knotfind_"+PDB;
         PrintWriter writer = new PrintWriter(trim_filename, "UTF-8"); 
         reader=new BufferedReader(new FileReader(file));
         tmp=new BufferedReader(new FileReader(file));
         String temp,temp1;
         while((temp=reader.readLine())!= null){
            String[] s=temp.split("\\s+");
             if(s[0].equals("ATOM")){
                if(s[2].equals("CA")){
                    if(s[4].equals("X")){
                        int residue = Integer.parseInt(s[5]);
                            if(residue == first_atom || list.contains(residue) || residue == last_atom){
                                writer.println(temp);
                            }
                        }
                    }
                }
            }
        
        writer.close();
        //System.out.println("Done Creating file at :" + PATH+PDB);
        //System.exit(0);
      }catch(IOException e){
      }
    }
}

class Res{
   double x;
   double y;
   double z;
   int index;
}

class Triangle{
   int res1;
   double distance;
}

class Plane{
   double a,b,c,d;
   int residual;
   public Plane(List<Res> res, int res1){
      Res p1=res.get(res1);
      Res p2=res.get(res1+1);
      Res p3=res.get(res1+2);
      if(res1<KnotFindOriginal.res.size()-2){
         this.residual=res1;
         a=(p1.y * (p2.z - p3.z)) + (p2.y) * (p3.z - p1.z) + (p3.y) * (p1.z - p2.z);

         b=(p1.z * (p2.x - p3.x)) + (p2.z) * (p3.x - p1.x) + (p3.z) * (p1.x - p2.x);

         c=(p1.x * (p2.y - p3.y)) + (p2.x) * (p3.y - p1.y) + (p3.x) * (p1.y - p2.y);

         d= - ((p1.x * ((p2.y * p3.z ) - (p3.y * p2.z )))
                 + (p2.x * ((p3.y * p1.z ) - (p1.y * p3.z )))
                 + (p3.x * ((p1.y * p2.z ) - (p2.y * p1.z ))));
      }
   }
}