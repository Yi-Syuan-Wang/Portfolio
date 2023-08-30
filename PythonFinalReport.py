import cv2
import math

joints = ['Head','T1','Hip','Knee','Ankle','Foot'] 
dots = []

def show_xy(event,x,y,flags,param): 
    if event == 1:
        dots.append([x, y])
        num = len(dots) 
        if num < 7:
            cv2.circle(img, (x, y), 10, (42,42,165), -1)
            cv2.putText(img, joints[num-1],(x+30,y+30), cv2.FONT_HERSHEY_SIMPLEX, 1, (42,42,165), 2, cv2.LINE_AA)
            if num > 1: 
                x1 = dots[num-2][0]
                y1 = dots[num-2][1]
                x2 = dots[num-1][0]
                y2 = dots[num-1][1]
                cv2.line(img,(x1,y1),(x2,y2),(139,236,255),2) 
                if num > 2:
                    x0 = dots[num-3][0]
                    y0 = dots[num-3][1]
                    print('%s: %d degrees'%((joints[num-2]),(angle(x0,y0,x1,y1,x2,y2))))
        cv2.imshow('Figure', img)
        cv2.imwrite('StickFigure.jpg', img)
        return

def angle(x0,y0,x1,y1,x2,y2):
    dx1 = x1-x0
    dy1 = y1-y0
    dx2 = x2-x1
    dy2 = y2-y1
    angle1 = math.atan2(dy1, dx1)
    angle1 = int(angle1 * 180/math.pi)
    angle2 = math.atan2(dy2, dx2)
    angle2 = int(angle2 * 180/math.pi)
    if angle1*angle2 >= 0:
        angleR = 180 - abs(angle1-angle2)
    else:
        angleR = abs(angle1) + abs(angle2)
        if angleR > 180:
            angleR = 360 - angleR
    return angleR


text ='''使用說明：
此程式用以繪製人體側面火柴圖。
繪製部位包含頭部(Head)、胸椎第一節(T1)、髖(Hip)、膝(Knee)、踝(Ankle)和足部(Foot)。'
並可返回名為「StickFigure」圖檔和各部位之相對夾角。
繪製完成後按任意鍵即可跳出視窗
注意：圖檔大小需500X500像素以上'''
print('*'*80)
print(text)
print('*'*80)


try:
    x = input('請輸入圖檔名稱(含副檔名): ')
    img = cv2.imread(x)
    #if img.shape[0] < 500 or img.shape[0] < 500:
        #print('圖檔需500X500像素以上')
    #else:
    cv2.imshow('Figure', img)
    cv2.setMouseCallback('Figure', show_xy)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
except AttributeError:
    print('例外發生，AttributeError')
except cv2.error:
    print('例外發生，cv2.error')
