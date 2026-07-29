package net.maizegenetics.stats.linearmodels;
/*
 * jGLiM: Java for General Linear Models
 * for more information: http://www.maizegenetics.net
 *
 * Copyright (C) 2005 Peter Bradbury
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

//package net.maizegenetics.jGLiM;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Dec 29, 2004
 * Time: 9:33:56 AM
 */
public interface Level extends Comparable<Level> {
    int getNumberOfSublevels();
    Comparable getSublevel(int sublevel);
    Comparable[] getSublevels();
    boolean contains(Comparable sublevel);

}
