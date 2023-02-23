logo <- sticker(imgurl,
        package = "Comparison\nof Microbial\nMethylated\nAdenines",
        #p_color = "#FC8EAC",
        #p_size = 24,
        #p_x = 1,
        #p_y = 1.4,
        url = "CoMMA",
        u_size = 24,
        u_color = "#FC8EAC",
        u_x = 1.15,
        u_y = 0.2,
        s_x = 1,
        s_y = 1,
        s_height = 0.5,
        h_color = "#FC8EAC",
        h_fill = "white",
        white_around_sticker = TRUE,
        dpi = 600,
        filename =
          "/Users/carlstone/Library/CloudStorage/Box-Box/Behringer_Lab_Box_Drive/Projects/LongTermExpEvo/MethylationProject/Rscripts/comma_logo.png")

logo <- logo + annotate("text", x = 1, y = 1, label = "CoMMA", color = "red")
