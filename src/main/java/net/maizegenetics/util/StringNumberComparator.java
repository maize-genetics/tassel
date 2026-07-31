/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and 
 * associations. For more information visit http://www.maizegenetics.net
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package net.maizegenetics.util;

import java.util.Comparator;
/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */

public class StringNumberComparator implements Comparator {

  public StringNumberComparator() {
  }

  public int compare(Object o1, Object o2) {
    try{
    Double d1=new Double(o1.toString());
    Double d2=new Double(o2.toString());
    return d1.compareTo(d2);
    }
    catch(NumberFormatException e)
      {String s1=o1.toString();
      String s2=o2.toString();
      return s1.compareTo(s2);
      }

  }
  /*
  public boolean equals(Object obj) {
    throw new java.lang.UnsupportedOperationException("Method equals() not yet implemented.");
  }
*/
}