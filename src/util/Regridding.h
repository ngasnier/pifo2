#pragma once

namespace pifo {
    class Regridding {
    private:
        // [ 0 1 2 ]      x
        // -------->
        // [ 2 1 0 ]
        static void optimizeGridIndices(double* x_in, long size_in, double* x_out, long size_out, bool cyclic, long* tab_i_in1, long* tab_i_in2, double* tab_x_adj1, double* tab_x_adj2);
        
    public:
        /**
         * Interpolation d'une grille 2D vers une autre par interpolation bilinéaire.
         * 
         * <p>Les valeurs à interpoler sont dans un tableau linéaire dont les 
         * coordonnées X et Y des axes sont fournies. Les intervalles séparant
         * les coordonnées X ou Y ne sont pas nécessairement réguliers.</p>
         * 
         * <p>Le tableau linéaire est ordonné de sorte que les indices pour un point 
         * de grille i, j sont calculés sous la forme k=j*NB coordonnéesY + X.</p>
         * 
         * <p>En sortie, on précise la liste des coordonnées X et Y que l'on
         * souhaite obtenir. Celles-ci peuvent être quelquonques, il doit y
         * en avoir autant que de point de sortie à obtenir.</p>
         * 
         * <p>Une interpolation bilinéaire entre les 4 points les plus proches
         * est effectuée pour obtenir les valeurs de sortie.</p>
         * 
         * <p>Si les coordonnées sont en dehors du domaine définit par les 
         * coordonnées de départ, la valeur du point le plus proche est prise.</p>
         * 
         * <p>Il est possible de faire reboucler la grille selon l'axe des X, 
         * dans ce cas toute coordonnée qui dépasse reboucle dans l'intervalle 
         * d'origine. L'intervalle dx entre le dernier point et le premier point 
         * est considéré comme équivalent à l'intervalle entre l'avant-dernier 
         * point et le dernier point.</p>
         * 
         * @param {type} x_in listing des coordonnées x d'origine. Doit être 
         * monotoniquement croissant ou décroissant.
         * @param {type} y_in listing des coordonnées y d'origine. Doit être 
         * monotoniquement croissant ou décroissant.
         * @param {type} data_in données d'entrée de taille x_in.length*y_in.length.
         * @param {type} cyclic est-ce que les données doivent être concidérées 
         * cycliques sur l'axe des X ? 
         * @param {type} x_out liste des coordonnées x de sorties souhaitées.
         * @param {type} y_out liste des coordonnées y de sorties souhaitées.
         * @param {type} data_out tableau de sortie, dans l'ordre des coordonnées
         * x_out, y_out fournies.
         * @returns {undefined}
         */
        static void bilinearRegrid(double* x_in, long in_width, 
                            double* y_in, long in_height,
                            double* data_in, long size_in, 
                            bool cyclic, 
                            double* x_out, double* y_out, 
                            double* data_out, long size_out);
    };
}