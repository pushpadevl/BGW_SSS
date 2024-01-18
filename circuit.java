import java.io.*;
import java.util.*;

// objective give the postfix expression
//given C = a + (b+c)*d-e
// enter circuit without spaces lowercase letters
// () > * > +
public class circuit {
    public static void main(String[] args) throws IOException{
        BufferedReader br  = new BufferedReader(new InputStreamReader(System.in));
        while(true){

            String str = br.readLine();
            Stack<Character> ST = new Stack<Character>();
            StringBuilder sb = new StringBuilder();
            //no spaces,  handle that later
            
            for(int i=0;i<str.length();i++){
                char cc = str.charAt(i);
                if(cc >= 'a' && cc <= 'z'){// handle more than 26 characters later
                    sb.append(cc);
                }else{
                    if(cc=='('){
                        ST.push(cc);
                    }else if(cc==')'){
                        while(ST.peek()!='(') sb.append(ST.pop()); 
                        ST.pop();
                        
                    }else if(cc=='+' ){
                        if(!ST.isEmpty() && ST.peek() == '*') sb.append(ST.pop());
                        ST.push(cc);
                    }else if(cc=='*'){
                        if(!ST.isEmpty() && ST.peek() == '*') sb.append(ST.pop());
                        ST.push(cc);
                    }
                }
                System.out.println(sb.toString());
            }
            
            if(!ST.isEmpty()){
                while(!ST.isEmpty()) sb.append(ST.pop());
            }
            System.out.println(sb.toString());
        }
    }    
}
